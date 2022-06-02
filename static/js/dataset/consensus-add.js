//
//	This test on the name change
//
$("#id-name-to-insert").on("change paste keyup", function () {
      $.ajax({
        url: $("#id-name-to-insert").attr("add-consensus-name-url"),
       // data: form.serialize(),
        // dataType: 'json',
        data : { 
        	consensus_name : $('#id-name-to-insert').val(),
			csrfmiddlewaretoken: '{{ csrf_token }}' }, // data sent with the post request
        success: function (data) {
          	$(document).find('#error_1_id_name').remove();
          	if (data.is_taken) {
             	var error_message = "<strong>" + data.error_message + "</strong>";
             	$("#id_dataset_name_error").append(
	            	$('<span/>', { 										// creates a dynamic div element on the fly
		        	id: 'error_1_id_name',
		        	class: 'col-10 fields_error',
		       		html: error_message,
	    		}));
				if ($('#id_consensus_fasta')[0].files.length > 0 && document.getElementById("id-name-to-insert").value.length > 0) {
					document.getElementById("id-save-button").disabled = true;
				}
				else{
					document.getElementById("id-save-button").disabled = false;
				}
         	}
			else{
				document.getElementById("id-save-button").disabled = false;
			}
        },
        
        // handle a non-successful response
        error : function(xhr,errmsg,err) {
            console.log(xhr.status + ": " + xhr.responseText); // provide a bit more info about the error to the console
        }
   });
});

$(document).ready(function(){
	// start with disable
	document.getElementById("id-save-button").disabled = true;
	
	// Prevent the form to be submitted on ENTER
	$('#id_form_consensus').submit(function(e){
	  	e.preventDefault();
	});
});

// open dialog upload consensus
$('#id_modal_add_new_consensus').on('click', function(){
	document.getElementById("id-save-button").disabled = true;
	document.getElementById("id-name-to-insert").value = '';
	document.getElementById("id_consensus_fasta").value = '';
});

// test file exists
function test_file_size() {
	if ($('#id_consensus_fasta')[0].files.length > 0 && document.getElementById("id-name-to-insert").value.length > 0) {
		document.getElementById("id-save-button").disabled = false;
    } else {
		document.getElementById("id-save-button").disabled = true;
    }
}

// add consensus
$(document).ready(function () {
	// alert($('#id_form_consensus').length)  // throws one if ID was found
	var form = $('#id_form_consensus')[0];
  	form.onsubmit = function(e){
	//$('#id_form_consensus').on('submit', function(e){
		e.preventDefault();
		
		//consensus name to insert
		$.ajax({
	        url: $('#id_form_consensus').attr("add-consensus-name-url"),
			type: "POST",
	        /// need to add crfs
	        data : { 
	        	consensus_name : $('#id-name-to-insert').val(),
				file_name: $('#id_consensus_fasta').serialize(),
	    		csrfmiddlewaretoken: '{{ csrf_token }}'
	        }, // data sent with the post request
	        		
	        success: function (data) {
	          if (data['is_ok']) {
	        		// add line with values
	        		$('#id_vaccine_status').append($('<option>', {
	        			value: data['value'],
	        			text: data['text'],
	        		 }));
	
		        	 /// add message with information
		        	 $('#id_messages_remove').append('<div class="alert alert-dismissible alert-success">' +
		        		data['message'] + '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
							'</div>');
					
					 /// Add new row in the table 
					 $('#id_tbody').append('<tr id="row_' + data['id'] + '" class="odd">' +
		                '<td class="name"><a href="#id_remove_modal" id="id_remove_dataset_modal" data-toggle="modal" title="Delete" ref_name="' +
						data['dataset_name'] + '" pk="' + data['id'] + '"><i class="fa fa-trash"></i> </a>' + data['dataset_name'] + '</td>' +
		                '<td class="last_change_date">Not set yet</td>' +
		                '<td class="creation_date">' + data['date_created'] + '</td>' +
		                '<td class="sequences"><span><i class="tip fa fa-info-circle" title="Consensus from projects: 0\n' +
		                'Consensus uploaded: 0\n' +
		                'References: 0"></i></span> (0/0/0) <a href="/datasets/datasets/' +
		                data['id'] + '/add_sequences_dataset" data-toggle="tooltip" title=""' +
		                'data-original-title="Add Consensus/References"><i class="fa fa-plus-square"></i> Add Seqs.</a></td>' +
		                '<td class="results"></td>' +
		                '</tr>')
	          }
	          else{
	        	  /// add message with information
	        	  $('#id_messages_remove').append('<div class="alert alert-dismissible alert-warning">' +
	        		data['message'] + '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
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
});


// remove dataset
$('#id-remove-button').on('click', function(){
	//id_name_to_insert
	$.ajax({
        url: $('#id-remove-button').attr("remove-dataset-name-url"),
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
        	  /// add message with information
        	  $('#id_messages_remove').append('<div class="alert alert-dismissible alert-success">' +
        		data['message'] + '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
				'</div>');
          }
          else{
        	  /// add message with information
        	  $('#id_messages_remove').append('<div class="alert alert-dismissible alert-warning">' +
        		data['message'] + '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
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

