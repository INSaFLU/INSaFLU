


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



