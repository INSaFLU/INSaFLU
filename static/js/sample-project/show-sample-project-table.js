//
//	This test on the name change
//

$('#collapseTwo_1').on('shown.bs.collapse', function () {
	show_variants_as_a_table();
});

//remove the tree from the screen
$('#collapseTwo_1').on('hidden.bs.collapse', function () {
	$('#showvariantsasatable').empty();
});

//draw phylocanvas
//set size of window "#phylocanvas" defined in 'static/css/flu-web-site.css'
// http://dmitrybaranovskiy.github.io/raphael/
function show_variants_as_a_table() {

 $.ajax({
 	/// spin 
 	beforeSend: function() {
 		$('#showvariantsasatable').empty();
 		$('#loader_showvariantsasatable').show();
 	},
 	complete: function(){
//  		$('#loader_showvariantsasatable').hide();
 	},
 	
	    data: { 
	    	project_id : $('#showvariantsasatable').attr("project_id"),
			csrfmiddlewaretoken: '{{ csrf_token }}'
	    }, // data sent with the get request
	    url: $('#showvariantsasatable').attr("show-variants-as-a-table-url"),
	    success: function (data) {
	    	if (data['is_ok']) {
	    	  
	    		/// create table
	    		d3.text(data['url_path_variant_table'], function(datasetText) {
	    			var parsedCSV = d3.tsv.parseRows(datasetText);

	    			var showvariantsasatable = document.getElementById('showvariantsasatable');
	    			    //.style("padding", "5px")
	    			    //.on("mouseover", function(){d3.select(this).style("background-color", "aliceblue")})
	    			    //.on("mouseout", function(){d3.select(this).style("background-color", "white")})
	    			    //.text(function(d){return d;});
	    			    //.style("font-size", "12px");
	    			
	    			var content = '<table id="table_with_variants_id">\n<thead>\n';
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
	    			$('#showvariantsasatable').append(content);
	    			
	    			var tf = new TableFilter('table_with_variants_id', {
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
	    			// ID	CHROM	POS	TYPE	REF	ALT	FREQ	COVERAGE	EVIDENCE	FTYPE	STRAND	NT_POS	AA_POS	EFFECT	NT CHANGE	AA CHANGE	AA CHANGE Alt	LOCUS_TAG	GENE	PRODUCT	VARIANTS IN INCOMPLETE LOCUS
	    			// covid_4	MN908947	514	snp	T	C	99.973	11115	C:11112 T:3	CDS	+	249/13203	83/4400	synonymous_variant	c.249T>C	p.His83His	p.H83H	orf1ab	orf1ab polyprotein	yes
	    			     // columns data types
	    			        col_types: [
	    			            'string',	// ID
	    			            'string',	// CHROM
	    			            'number',	// POS
	    			            'string',	// TYPE
	    			            'string',	// REF
	    			            'string',	// ALT
	    			            { type: 'formatted-number', decimal: ',', thousands: '.' }, // FREQ
	    			            //'formatted-number', // defaults to '.' for decimal and ',' for thousands
	    			            'number',	// COVERAGE
	    			            'string',	// EVIDENCE
	    			            'string',	// FTYPE
	    			            'string',	// STRAND  +/-
	    			            'string',	// NT_POS
	    			            'string',	// AA_POS
	    			            'string',	// EFFECT
	    			            'string',	// NT CHANGE
	    			            'string',	// AA CHANGE Alt
	    			            'string',	// AA CHANGE
	    			            'string',	// LOCUS_TAG
	    			            'string',	// GENE
	    			            'string',	// PRODUCT
	    			            'string',	// VARIANTS IN INCOMPLETE LOCUS [yes/empty]
	    			        ],
	    			        col_widths: [
	    			        	"7%",	// ID
	    			            "6%",	// CHROM
	    			            "3%",	// POS
	    			            "3%",	// TYPE
	    			            "3%",	// REF
	    			            "3%",	// ALT
	    			            "3%", // FREQ
	    			            "4%",	// COVERAGE
	    			            null,	// EVIDENCE
	    			            "3.5%",	// FTYPE
	    			            "3%",	// STRAND  +/-
	    			            null,	// NT_POS
	    			            null,	// AA_POS
	    			            null,	// EFFECT
	    			            "4%",	// NT CHANGE
	    			            null,	// AA CHANGE
	    			            null,	// AA CHANGE Alt
	    			            null,	// LOCUS_TAG
	    			            "3%",	// GENE
	    			            null,	// PRODUCT
	    			            null,	// VARIANTS IN INCOMPLETE LOCUS [yes/empty]
	    			        ],
	    			        
	    			        col_1: 'select',
	    			        col_3: 'select',
	    			        col_9: 'select',
	    			        col_10: 'select',
	    			        col_13: 'select',
	    			        col_17: 'select',
	    			        col_19: 'select',
	    			        // Sort extension: in this example the column data types are provided by the
	    			        // 'col_types' property. The sort extension also has a 'types' property
	    			        // defining the columns data type for column sorting. If the 'types'
	    			        // property is not defined, the sorting extension will fallback to
	    			        // the 'col_types' definitions.
	    			        extensions: [{ name: 'sort' }],
	    				});
	    				tf.init();

	    			// add download possibility
	    			content = '<a href="#" class="helpBtn" onclick="download_table_as_csv(\'table_with_variants_id\');">Download as CSV</a>'
	    			var element_rdiv = document.querySelector(".rdiv");
	    			if (typeof(element_rdiv) != "undefined"){
		    			let span = document.createElement('span');
		    			span.innerHTML = content;
		    			element_rdiv.appendChild(span);
	    			}

	    			// hide loading icon
	    			$('#loader_showvariantsasatable').hide();
	    		});
	      }
	      else{
	    	  $('#loader_showvariantsasatable').hide();
	    	  $('#showvariantsasatable').append('<div class="alert alert-warning alert-dismissable text-center"><strong>Fail</strong> to load the variants file.</div>')
	      }
	    },
	    
	    // handle a non-successful response
	    error : function(xhr,errmsg,err) {
	    	$('#loader_showvariantsasatable').hide();
	    	$('#showvariantsasatable').append('<div class="alert alert-warning alert-dismissable text-center"><strong>Fail</strong> to load the variants file.</div>')
	        console.log(xhr.status + ": " + xhr.responseText); // provide a bit more info about the error to the console
	    }
 });
}


//Quick and simple export target #table_id into a csv
// <a href="#" onclick="download_table_as_csv('my_id_table_to_export');">Download as CSV</a>
// <img src="http://127.0.0.1:8000/static/vendor/tablefilter/style/themes/blank.png" class="sort-arrow">
function download_table_as_csv(table_id, separator = ',') {
	
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

