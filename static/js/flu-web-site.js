


$().ready(function(){
	$( "a" ).hover(
	  function() {   
		   var title = $(this).attr("data-title");  			// extracts the title using the data-title attr applied to the 'a' tag
		   if (title != null){
			    $('<div/>', { 										// creates a dynamic div element on the fly
			        text: title,
			        class: 'box_tooltip_flu'
			    }).appendTo(this);  								// append to 'a' element
		     }
		},
		function() {
		    $(document).find("div.box_tooltip_flu").remove(); 	// on hover out, finds the dynamic element and removes it.
		}
	);
});


$().ready(function(){
    $('#my_modal_runner').DjangoModalRunner({
    	 on_show_modal: function(){
            alert('show');
        },
        on_hide_modal: function(){
        	alert('hide without action');
//        	$('#my_modal_runner').append(
//	        	$('<div/>', { 										// creates a dynamic div element on the fly
//			        text: 'vdsddsvdstitle',
//			        class: 'box_tooltip_flu_'
//		    	})) 
        },
        on_hide_modal_after_submit: function(){
            alert('hide with action');
        },
        on_submit: function(){
            alert('submit');
        },
        on_done: function(){
            alert('done');
        }
    });
});

/// Call modal for VaccineStatus AddSample
$().ready(function(){
    $('#vaccine_add_modal').DjangoModalRunner({
    	 on_show_modal: function(){
            alert('show');
        },
        on_hide_modal: function(){
        	alert('hide without action');
        $('#id_vaccine_status').append(
	        	$('<option/>', { 										// creates a dynamic div element on the fly
			        value: 'vdsddsvdstitle',
			        text: 'vdsddsvdstitle',
		    	}))
        },
        on_hide_modal_after_submit: function(){
            alert('hide with action');
        },
        on_submit: function(){
            alert('submit');
        },
        on_done: function(){
            alert('done');
        }
    });
});

/// Call modal for DatSet AddSample
$().ready(function(){
    $('#data_set_add_modal').DjangoModalRunner({
    	 on_show_modal: function(){
            alert('show');
        },
        on_hide_modal: function(){
        	alert('hide without action');
        
        	$('#id_data_set').append(
	        	$('<option/>', { 										// creates a dynamic div element on the fly
			        value: 'vdsddsvdstitle',
			        text: 'vdsddsvdstitle',
		    	}))
		    
            
        },
        on_hide_modal_after_submit: function(){
            alert('hide with action');
        },
        on_submit: function(){
            alert('submit');
        },
        on_done: function(){
            alert('done');
        }
    });
});



/// Add class 'rounded-group-box-legend' on the fly
$().ready(function(){
	$("legend").addClass('rounded-group-box-legend');
});

/// change help-block class by fields_error
$(document).ready(function(){
    $('.help-block').removeClass('help-block').addClass('fields_error');
});

/// made a click in the nav tab, this is a caveat
$(document).ready(function(){
    $('.nav li.active').removeClass('active').find('a').trigger('click');
});

/// add function to toggle the checkbox in tables
document.getElementById("checkBoxAll").addEventListener ("click", toggle_check_box_all, false);
function toggle_check_box_all(source) {
	var remember = document.getElementById('checkBoxAll');
    checkboxes = document.getElementsByName('select_ref');
    for(var i=0, n=checkboxes.length;i<n;i++) {
		checkboxes[i].checked = remember.checked;
	}
	
	$.ajax({
		url: $('#table_with_check_id').attr("set-check-box-values-url"),
		data : { check_box_all : remember.checked }, // data sent with the post request
		success: function (data) { },
	});
};

$(document).ready(function(){
	var elements = document.getElementsByName("select_ref");
	for(var i = 0, n = elements.length; i < n; i++){
		elements[i].addEventListener('click', toggle_check_box, false);
	}
	// set all checked in the server
	$.ajax({
		url: $('#table_with_check_id').attr("set-check-box-values-url"),
		data : { get_check_box_single : '1' }, // data sent with the post request
		success: function (data) {
			for (key in data){
				if (key === 'is_ok'){ continue; }
				var remember = document.getElementById(key);
				if (remember != null){
					remember.checked = data[key];
				}
			}
		},
	});
});

function toggle_check_box(source) {
	$.ajax({
		url: $('#table_with_check_id').attr("set-check-box-values-url"),
		data : { check_box : source.srcElement.checked,
				 value : source.srcElement.value }, // data sent with the post request
		success: function (data) { },
	});
};


//// everything about checkBox
$(document).ready(function(){
	var remember = document.getElementById('checkBoxAll');
    var check_box_all_session = $('#table_with_check_id').attr("check_box_all");
    if (check_box_all_session == "true" ){
    	remember.checked = true;
    	checkboxes = document.getElementsByName('select_ref');
    	for(var i=0, n=checkboxes.length;i<n;i++) {
    		checkboxes[i].checked = remember.checked;
		}
	}
	else{
		remember.checked = false;
	}
});

// if the user pressed the 
$(function() { //shorthand document.ready function
	$('#id_add_all_checked').on('submit', function (e) {
		e.preventDefault();  //prevent form from submitting
	    $.ajax({
			url: $('#table_with_check_id').attr("set-check-box-values-url"),
	        data : { count_check_boxes : '1' }, // data sent with the post request
			success: function (data) {
				if (data['count_check_boxes'] < 1){
					alert('There is no samples selected.\nPlease, for this option select some samples.');
					return false;
				}
			},
			// handle a non-successful response
	        error : function(xhr,errmsg,err) {
	            alert(errmsg);
	            console.log(xhr.status + ": " + xhr.responseText); // provide a bit more info about the error to the console
	        },
		});
		return true;
	});
});  	



