//
//	This test on the name change
//

$('#collapseAln2Pheno').on('shown.bs.collapse', function () {
	show_aln2pheno();
});

//remove the tree from the screen
$('#collapseAln2Pheno').on('hidden.bs.collapse', function () {
	$('#showaln2pheno').empty();
});

//draw phylocanvas
//set size of window "#phylocanvas" defined in 'static/css/flu-web-site.css'
// http://dmitrybaranovskiy.github.io/raphael/
function show_aln2pheno() {

 $.ajax({
 	/// spin 
 	beforeSend: function() {
 		$('#showaln2pheno').empty();
 		$('#loader_showaln2pheno').show();
 	},
 	complete: function(){
//  		$('#loader_showvariantsasatable').hide();
 	},
 	
	    data: { 
	    	project_id : $('#showaln2pheno').attr("project_id"),
			csrfmiddlewaretoken: '{{ csrf_token }}'
	    }, // data sent with the get request
	    url: $('#showaln2pheno').attr("show-aln2pheno-url"),
	    success: function (data) {
	    	if (data['is_ok']) {
	    	  
	    		/// create table
	    		d3.text(data['url_path_aln2pheno'], function(datasetText) {
	    			var parsedCSV = d3.tsv.parseRows(datasetText);

	    			var showaln2pheno = document.getElementById('showaln2pheno');
	    			
	    			var content = '<table id="table_with_aln2pheno_id">\n<thead>\n';
	    			var count = 0;
	    			var number_cells = 0;
	    			var number_cells_processed = 0;
	    			parsedCSV.forEach(function(row) {
	    				content += "<tr>";
	    				number_cells_processed = 0;
		    	        row.forEach(function(cell) {
		    	            content += "<td>" + cell + "</td>" ;
		    	            number_cells_processed += 1
		    	            if (count === 0){
		    	            	number_cells += 1
		    	            }
		    	        });
		    	        /// sometimes the lines doesn't have all fields
		    	        if (count != 0 && number_cells_processed != number_cells){
		    	        	for (i = 0; i < (number_cells - number_cells_processed); i++) {
		    	        		content += "<td></td>";
		    	        	}
		    	        }
		    	        content += "</tr>";
		    	        
	    				if (count === 0){
	    					content += "</thead>\n<tbody>"
	    				}
	    				count += 1;
		    	    });
	    			content += '</tbody></table>'
//	    			console.log(content);
	    			$('#showaln2pheno').append(content);
	    			
	    			var tf = new TableFilter('table_with_aln2pheno_id', {
	    					base_path: data['static_table_filter'] + '/',
	    					auto_filter: {
	    			            delay: 500 //milliseconds
	    			        },
	    					alternate_rows: true,
	    			        rows_counter: true,
	    			        btn_reset: true,
	    			   //     status_bar: true,
	    			        msg_filter: 'Filtering...',
	    			        rows_counter: {
	    			            text: 'Variants: '
	    			        },
	    					grid_layout: {
	    						width: '100%'
	    					},
	    			// Sequence	NTD	RBD	RBM	RBD_RBM	S2	FP	S1	S1_other	SP	Flagged mutations	All mutations	Nflagged	Nmutations	lineage
	    			     // columns data types
	    			        col_types: [
	    			            'string',	// Sequence
	    			            'string',	// lineage
	    			            'string',	// NTD
	    			            'string',	// RDB
	    			            'string',	// S2
	    			            'string',	// F1
	    			            'string',	// S1
	    			            'string',	// S1_other
	    			            'string',	// SP
	    			            'string',	// Flagged mutations
	    			            'string',	// All mutations
	    			            'number',	// Nflagged
	    			            'number',	// Nmutations
	    			        ],
	    			        col_widths: [
	    			        	"10%",	// Sequence
								"7%",	// lineage								
	    			            "7%",	// NTD
	    			            "7%",	// RDB
	    			            "7%",	// S2
	    			            "7%",	// F1
	    			            "7%",   // S1
	    			            "7%",	// S1_other
	    			            "7%",	// SP
	    			            "7%",	// Flagged mutations
	    			            null,	// All mutations
	    			            "7%",	// Nflagged
	    			            "7%",	// Nmutations
	    			        ],
	    			        
	    			        col_1: 'select',
	    			        col_2: 'select',
	    			        col_3: 'select',
	    			        col_4: 'select',
	    			        col_5: 'select',
	    			        col_6: 'select',
	    			        col_7: 'select',
                            col_8: 'select',
                            col_9: 'select',
	    			        // Sort extension: in this example the column data types are provided by the
	    			        // 'col_types' property. The sort extension also has a 'types' property
	    			        // defining the columns data type for column sorting. If the 'types'
	    			        // property is not defined, the sorting extension will fallback to
	    			        // the 'col_types' definitions.
	    			        extensions: [{ name: 'sort' }],
	    				});
	    				tf.init();

	    			// add download possibility
	    			content = '<a href="#" class="helpBtn" onclick="download_aln2pheno_as_csv(\'table_with_aln2pheno_id\');">Download as CSV</a>'
	    			var element_rdiv = document.querySelector(".rdiv");
	    			if (typeof(element_rdiv) != "undefined"){
		    			let span = document.createElement('span');
		    			span.innerHTML = content;
		    			element_rdiv.appendChild(span);
	    			}

	    			// hide loading icon
	    			$('#loader_showaln2pheno').hide();
	    		});
	      }
	      else{
	    	  $('#loader_showaln2pheno').hide();
	    	  $('#showaln2pheno').append('<div class="alert alert-warning alert-dismissable text-center"><strong>Fail</strong> to load the variants file.</div>')
	      }
	    },
	    
	    // handle a non-successful response
	    error : function(xhr,errmsg,err) {
	    	$('#loader_showaln2pheno').hide();
	    	$('#showaln2pheno').append('<div class="alert alert-warning alert-dismissable text-center"><strong>Fail</strong> to load the variants file.</div>')
	        console.log(xhr.status + ": " + xhr.responseText); // provide a bit more info about the error to the console
	    }
 });
}


//Quick and simple export target #table_id into a csv
// <a href="#" onclick="download_table_as_csv('my_id_table_to_export');">Download as CSV</a>
// <img src="http://127.0.0.1:8000/static/vendor/tablefilter/style/themes/blank.png" class="sort-arrow">
function download_aln2pheno_as_csv(table_id, separator = ',') {
	
	/// do the Head, Select rows from table_id
	var rows = document.querySelectorAll('div.grd_headTblCont tr');
	// Construct csv, 
	var csv = [];
	for (var i = 0; i < rows.length; i++) {
        var row = [], cols = rows[i].querySelectorAll('td, th');
        for (var j = 0; j < cols.length; j++) {
        	 var data = cols[j].innerText.replace(/(\r\n|\n|\r)/gm, '').replace(/(\s\s)/gm, ' ');
        	 data = data.replace(/"/g, '""');
        	 row.push('"' + data + '"');
        	 break;
        }
        csv.push(row.join(separator));
    }
    
    /// do the body, Select rows from table_id
    rows = document.querySelectorAll('table#' + table_id + ' tr.odd, table#' + table_id + ' tr.even');
    for (var i = 0; i < rows.length; i++) {
        var row = [], cols = rows[i].querySelectorAll('td, td');
        for (var j = 0; j < cols.length; j++) {
            // Clean innertext to remove multiple spaces and jumpline (break csv)
            var data = cols[j].innerText.replace(/(\r\n|\n|\r)/gm, '').replace(/(\s\s)/gm, ' ');
            // Escape double-quote with double-double-quote (see https://stackoverflow.com/questions/17808511/properly-escape-a-double-quote-in-csv)
            data = data.replace(/"/g, '""');
            // Push escaped string
            row.push('"' + data + '"');
        }
        csv.push(row.join(separator));
    }
    var csv_string = csv.join('\n');
    // Download it
    var filename = 'export_' + table_id + '_' + new Date().toLocaleDateString() + '.csv';
    var link = document.createElement('a');
    link.style.display = 'none';
    link.setAttribute('target', '_blank');
    link.setAttribute('href', 'data:text/csv;charset=utf-8,' + encodeURIComponent(csv_string));
    link.setAttribute('download', filename);
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
}

