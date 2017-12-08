


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



