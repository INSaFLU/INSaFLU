//
//	This test on the name change
//
//$("#phylocanvas").on("change paste keyup", function () {

$('#collapseSixth').on('shown.bs.collapse', function () {
	show_coverage_table();
});

//remove the tree from the screen
$('#collapseSixth').on('hidden.bs.collapse', function () {
	$('#show_coverage_id').empty();
});

//draw phylocanvas
//set size of window "#phylocanvas" defined in 'static/css/flu-web-site.css'
function show_coverage_table() {

 $.ajax({
 	/// spin 
 	beforeSend: function() {
 		$('#show_coverage_id').empty();
 		$('#loader_coverage_id').show();
 	},
 	complete: function(){
//  		$('#loader_coverage_id').hide();
 	},
 	
	    data: { 
	    	project_id : $('#show_coverage_id').attr("project_id"),
	    	client_width : document.getElementById('headingSixth').offsetWidth,
			csrfmiddlewaretoken: '{{ csrf_token }}'
	    }, // data sent with the get request
	    url: $('#show_coverage_id').attr("show-coverage-as-a-table-url"),
	    success: function (data) {
	    	if (data['is_ok']) {
	    	  
	    		// add content
	    		$('#show_coverage_id').append(data['content']);

	    		// get size of table
	    		var number_coluns = $("#table_with_coverage_variants_id").attr("number_coluns");
	    		var number_rows = $("#table_with_coverage_variants_id").attr("number_rows");
	    		
	    		//// draw bars
	    		var y_size = 25;
	    		for (i = 1; i <= number_rows; i++) {
	    			for (x = 1; x <= number_coluns; x++) {
	    				var id_to_search = 'id_table-coverage_' + i + '_' + x;
	    				var id_to_get_value = 'id_table-coverage_content_' + i + '_' + x;
	    				var value_coverage = $("#" + id_to_get_value).attr("value_data");
	    				var value_limit_coverage = $("#" + id_to_get_value).attr("value_limit_coverage");
	    				var value_data_average = $("#" + id_to_get_value).attr("value_data_average");
	    				var color_graphic = $("#" + id_to_get_value).attr("color_graphic");
	    				var x_size = $("#" + id_to_get_value).attr("size");
	    				var draw_x = x_size * value_coverage / 100;
	    				
	    				var paper = Raphael(id_to_search, x_size, y_size);
	    				if (draw_x > 1){
				            var rectangle = paper.rect(0, 0, draw_x, y_size, 5);
				            var filler = { fill: color_graphic, cursor: 'pointer' };
				            rectangle.attr(filler);
	    				}
	    				paper.text(x_size / 2 - 3, 12, value_data_average + ' - ' + value_coverage + '%' + 
	    					' (' + value_limit_coverage + 'x)').attr({
			            	fill: '#000000',
			            	"font-family": "Arial",
			            	"font-size":"15px",
			            	"font-weight": "normal",
			            });
	    			}
	    		}
	            
    			// hide loading icon
    			$('#loader_coverage_id').hide();
	      }
	      else{
	    	  $('#loader_coverage_id').hide();
	    	  $('#show_coverage_id').append('<div class="alert alert-warning alert-dismissable text-center"><strong>Fail</strong> to load coverages.</div>')
	      }
	    },
	    
	    // handle a non-successful response
	    error : function(xhr,errmsg,err) {
	    	$('#loader_coverage_id').hide();
	    	$('#show_coverage_id').append('<div class="alert alert-warning alert-dismissable text-center"><strong>Fail</strong> to load coverages.</div>')
	        console.log(xhr.status + ": " + xhr.responseText); // provide a bit more info about the error to the console
	    }
 });
}

