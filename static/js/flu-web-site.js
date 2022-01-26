


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
/*$(document).ready(function(){
    $('.nav li.active').removeClass('active').find('a').trigger('click');
});*/

/// set wait screen
function wait_screen() {
  $.blockUI({ css: { 
	    border: 'none', 
	    padding: '15px', 
	    backgroundColor: '#000', 
	    '-webkit-border-radius': '10px', 
	    '-moz-border-radius': '10px', 
	    opacity: .5, 
	    color: '#fff' 
	} });
}
//	START expanding logo
//window.onscroll = function() {
//	growShrinkLogo()
//};

// https://stackoverflow.com/questions/24765155/shrinking-navigation-bar-when-scrolling-down-bootstrap3
//function growShrinkLogo() {
//	var Logo = document.getElementById("main_logo_id_")
//	if (document.body.scrollTop > 5 || document.documentElement.scrollTop > 5) {
//		Logo.style.width = '30px';
//	} else {
//		Logo.style.width = '60px';
//	}
//}
//	END expanding logo