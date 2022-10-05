//
//	This test on the name change
//
$("#id-name-to-insert").on("change paste keyup", function () {
      $.ajax({
        url: $("#id-name-to-insert").attr("add-dataset-name-url"),
       // data: form.serialize(),
        // dataType: 'json',
        data : { 
        	dataset_name : $('#id-name-to-insert').val(),
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
				document.getElementById("id-save-button").disabled = true;
			 
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

// add sample
$('#id-save-button').on('click', function(){
	//id_name_to_insert
	$.ajax({
        url: $('#id-save-button').attr("add-dataset-name-url"),
        /// need to add crfs
        data : { 
        	dataset_name : $('#id-name-to-insert').val(),
    		csrfmiddlewaretoken: '{{ csrf_token }}'
        }, // data sent with the post request
        		
        success: function (data) {
          if (data['is_ok']) {

				// remove empty text in table if exists
				var element = document.getElementById('id_empty_text_table')
				if (element != null){
					element.remove()
				}
			
	        	 /// add message with information
	        	 $('#id_messages_remove').append('<div class="alert alert-dismissible alert-success">' +
	        		data['message'] + '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
						'</div>');
				
				 /// Add new row in the table 
				 $('#id_tbody').append('<tr id="row_' + data['id'] + '" class="odd">' +
	                '<td class="name"><a href="#id_remove_modal" id="id_remove_dataset_modal" data-toggle="modal" title="Delete" ref_name="' +
					data['dataset_name'] + '" pk="' + data['id'] + '"><i class="fa fa-trash"></i> </a>' + data['dataset_name'] + '</td>' +
					'<td class="build">Generic</td>' +
	                '<td class="last_change_date">Not set yet</td>' +
	                '<td class="creation_date">' + data['date_created'] + '</td>' +
					'<td class="totla_alerts">0</td>' +
	                '<td class="sequences"><span><i class="tip fa fa-info-circle" title="Consensus from projects: 0\n' +
	                'Consensus uploaded: 0\nReferences: 0"></i></span> (0/0/0)' +
	                '<div class="btn-group"><button type="button" class="btn btn-primary dropdown-toggle" data-toggle="dropdown" aria-haspopup="true" ' +
					'aria-expanded="false" title="Add Consensus/References/Projects" style="margin-left: 8px;">Add Sequences</button>' +
					'<div class="dropdown-menu" x-placement="bottom-start" style="position: absolute; ' +
					'transform: translate3d(0px, 38px, 0px); top: 0px; left: 0px; will-change: transform;">' +
					'<a rel="nofollow" class="dropdown-item" href="/datasets/datasets/' + data['id'] + '/add_references_dataset"> ' +
					'Add References</a><a rel="nofollow" class="dropdown-item" href="/datasets/datasets/' + data['id'] + '/add_projects_dataset">' +
					' Add Consensus from Projects</a><a rel="nofollow" class="dropdown-item" href="/datasets/datasets/' + data['id'] + '/add_consensus_dataset">' +
					' Add your own Consensus</a></div>' +
					'<td class="results"><a rel="nofollow" href="/datasets/' + data['id'] + '/show_dataset_settings" data-toggle="tooltip" title="Dataset settings">' +
					'<span ><i class="fa fa-magic padding-button-table"></i></span></a></td>' +
					'</div>' +
	                '</td></tr>')
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

