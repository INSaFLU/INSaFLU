//
//	This test on the name change
//
//$("#phylocanvas").on("change paste keyup", function () {

$('#collapseTwo').on('shown.bs.collapse', function () {
	draw_phylo_canvas();
});

//remove the tree from the screen
$('#collapseTwo').on('hidden.bs.collapse', function () {
	$('#phylocanvas').empty();
});

//catch the element from tree nucleotides combo box
$("#combo_select_elements_phylocanvas_id").change(function () {
	draw_phylo_canvas();
});


// draw phylocanvas
function draw_phylo_canvas() {
	var element_selected = $('#combo_select_elements_phylocanvas_id option:selected').val();
    
    $.ajax({
    	/// spin 
    	beforeSend: function() {
    		$('#phylocanvas').empty();
    		$('#loader_phylocanvas').show();
    	},
    	complete: function(){
    		$('#loader_phylocanvas').hide();
    	},
    	
	    data : { 
	    	project_id : $('#phylocanvas').attr("project_id"),
	    	key_element_name : element_selected,
	    }, // data sent with the get request
	    
	    url: $('#phylocanvas').attr("show-phylo-canvas-url"),
	    success: function (data) {
	      if (data['is_ok']) {
	    	  (function (Phylocanvas) {
	    	      var tree = Phylocanvas.createTree('phylocanvas', { history: false, });
	    	      tree.load(data['tree']);
	    	      tree.setTreeType('rectangular');
	    	  })(window.Phylocanvas);
	    	  
	    	  //// set the files
			  $('#tree_nwk_id').empty()
			  $('#tree_nwk_id').append(data['tree_nwk_id'])
			  $('#tree_tree_id').empty()
			  $('#tree_tree_id').append(data['tree_tree_id'])
	      }
	      else{
	    	  $('#phylocanvas').append('<div class="alert alert-warning alert-dismissable"><strong>Fail</strong> to load the tree.</div>')
	      }
	    },
	    
	    // handle a non-successful response
	    error : function(xhr,errmsg,err) {
	    	$('#phylocanvas').append('<div class="alert alert-warning alert-dismissable"><strong>Fail</strong> to load the tree.</div>')
	        console.log(xhr.status + ": " + xhr.responseText); // provide a bit more info about the error to the console
	    }
    });
}

/**
 * SHOW nucleotides alignments 
 *
 */
$('#collapseThree').on('shown.bs.collapse', function () {
	draw_nucleotide_alignments();
});

//remove the tree from the screen
$('#collapseThree').on('hidden.bs.collapse', function () {
	$('#msa_viewer_nucleote_id').empty();
});

//catch the element from tree nucleotides combo box
$("#combo_select_elements_nucleotids_alignments_id").change(function () {
	draw_nucleotide_alignments();
});

//draw nucleotide alignments
function draw_nucleotide_alignments() {
	var element_selected = $('#combo_select_elements_nucleotids_alignments_id option:selected').val();
    
    $.ajax({
    	/// spin 
    	beforeSend: function() {
    		$('#msa_viewer_nucleote_id').empty();
    		$('#loader_msa_viewer_nucleote_id').show();
    	},
    	complete: function(){
    		$('#loader_msa_viewer_nucleote_id').hide();
    	},
    	
	    data : { 
	    	project_id : $('#msa_viewer_nucleote_id').attr("project_id"),
	    	key_element_name : element_selected,
	    }, // data sent with the get request
	    
	    url: $('#msa_viewer_nucleote_id').attr("show-msa-nucleotide-url"),
	    success: function (data) {
	      if (data['is_ok']) {
	    	  
	    	  data['is_ok']
	    	  var rootDiv = document.getElementById("msa_viewer_nucleote_id");
	    	  var labelNameLength = data['max_length_label'] * 13 + 10;	// multiply by 13 because the labelFontsize is 13; https://github.com/wilzbach/msa in zoomer section
	    	  var m = msa({
	    			el: rootDiv,
	    			vis: {
	    			      labelId: false
	    			},
	    	        zoomer: {
	    	        	labelNameLength: labelNameLength,
	    	        },
	    	  });
	    	  msa.io.fasta.read(data['alignment_fasta_show_id'], function(err, seqs){
	    			m.seqs.reset(seqs);
	    			m.render();
	    	  });
	    	  
	    	  //// set the files
        	  $('#alignment_fasta_id').empty()
			  $('#alignment_fasta_id').append(data['alignment_fasta_id'])
			  $('#alignment_nex_id').empty()
			  $('#alignment_nex_id').append(data['alignment_nex_id'])
	      }
	      else{
	    	  $('#msa_viewer_nucleote_id').append('<div class="alert alert-warning alert-dismissable"><strong>Fail</strong> to load the alignment.</div>')
	      }
	    },
	    
	    // handle a non-successful response
	    error : function(xhr,errmsg,err) {
	    	$('#msa_viewer_nucleote_id').append('<div class="alert alert-warning alert-dismissable"><strong>Fail</strong> to load the alignment.</div>')
	        console.log(xhr.status + ": " + xhr.responseText); // provide a bit more info about the error to the console
	    }
    });
}

/**
 * SHOW protein alignments 
 *
 */
$('#collapseFourth').on('shown.bs.collapse', function () {
	draw_protein_alignments();
});

//remove the tree from the screen
$('#collapseFourth').on('hidden.bs.collapse', function () {
	$('#msa_viewer_amino_id').empty();
});

//catch the element from tree nucleotides combo box
$("#combo_select_elements_amino_alignments_id").change(function () {
	// change elements in combo_select_gene_amino_alignments_id
	var element_selected = $('#combo_select_elements_amino_alignments_id option:selected').val();
	$.ajax({
    	/// spin 
    	beforeSend: function() {
    		$('#combo_select_gene_amino_alignments_id').empty();
    	},
    	complete: function(){
    		draw_protein_alignments();
    	},
    	data : { 
	    	project_id : $('#msa_viewer_amino_id').attr("project_id"),
	    	key_element_name : element_selected,
	    },
    	url: $('#combo_select_elements_amino_alignments_id').attr("get-cds-from-element-url"),
    	success: function (data) {
  	      if (data['is_ok']) {
  	    	  for( i = 0; i < data['vect_genes'].length; i++ ){ 
	  	    	  if (i == 0){
	  	    		  $('#combo_select_gene_amino_alignments_id').append('<option selected value="' + data['vect_genes'][i] + '">' + data['vect_genes'][i] + '</option>');
	  	    	  }
	  	    	  else{
	  	    		  $('#combo_select_gene_amino_alignments_id').append('<option value="' + data['vect_genes'][i] + '">' + data['vect_genes'][i] + '</option>');
	  	    	  }
  	    	  }
  	      }
    	}
	});
	
});

//catch the element from tree nucleotides combo box
$("#combo_select_gene_amino_alignments_id").change(function () {
	draw_protein_alignments();
});


//draw nucleotide alignments
function draw_protein_alignments() {
	var element_selected = $('#combo_select_elements_amino_alignments_id option:selected').val();
	var gene_selected = $('#combo_select_gene_amino_alignments_id option:selected').val();
    
    $.ajax({
    	/// spin 
    	beforeSend: function() {
    		$('#msa_viewer_amino_id').empty();
    		$('#loader_msa_viewer_amino_id').show();
    	},
    	complete: function(){
    		$('#loader_msa_viewer_amino_id').hide();
    	},
    	
	    data : { 
	    	project_id : $('#msa_viewer_amino_id').attr("project_id"),
	    	key_element_name : element_selected,
	    	key_gene_name : gene_selected,
	    }, // data sent with the get request
	    
	    url: $('#msa_viewer_amino_id').attr("show-msa-amino-url"),
	    success: function (data) {
	      if (data['is_ok']) {
	    	  
	    	  var rootDiv = document.getElementById("msa_viewer_amino_id");
	    	  var labelNameLength = data['max_length_label'] * 14 + 10;	// multiply by 14 because the labelFontsize is 13; https://github.com/wilzbach/msa in zoomer section
	    	  var m = msa({
	    			el: rootDiv,
	    			vis: {
	    			      labelId: false
	    			},
	    	        zoomer: {
	    	        	labelNameLength: labelNameLength,
	    	        },
	    	  });
	    	  msa.io.fasta.read(data['alignment_amino_fasta_show_id'], function(err, seqs){
	    			m.seqs.reset(seqs);
	    			m.render();
	    	  });
	    	  
	    	  //// set the files
        	  $('#alignment_amino_fasta_id').empty()
			  $('#alignment_amino_fasta_id').append(data['alignment_amino_fasta_id'])
			  $('#alignment_amino_nex_id').empty()
			  $('#alignment_amino_nex_id').append(data['alignment_amino_nex_id'])
	      }
	      else{
	    	  $('#msa_viewer_amino_id').append('<div class="alert alert-warning alert-dismissable"><strong>Fail</strong> to load the alignment.</div>')
	      }
	    },
	    
	    // handle a non-successful response
	    error : function(xhr,errmsg,err) {
	    	$('#msa_viewer_amino_id').append('<div class="alert alert-warning alert-dismissable"><strong>Fail</strong> to load the alignment.</div>')
	        console.log(xhr.status + ": " + xhr.responseText); // provide a bit more info about the error to the console
	    }
    });
}

// graphics to show count variations
$('#collapseFifth').on('shown.bs.collapse', function () {
	if (! $(".chartjs-size-monitor")[0]){
		draw_count_variations_chart();
	}
});

//remove the tree from the screen
$('#collapseFifth').on('hidden.bs.collapse', function () {
	$('#show_count_variations_id').empty();
});


/// get count variations from server and show in a graph
function draw_count_variations_chart(){
	
	
	$.ajax({
    	/// spin 
    	beforeSend: function() {
    		$('#show_count_variations_id').empty();
    		$('#loader_count_variations_id').show();
    	},
    	complete: function(){
    		$('#loader_count_variations_id').hide();
    	},
    	
	    data : { 
	    	project_id : $('#canvas_count_variations_id').attr("project_id"),
	    }, // data sent with the get request
	    
	    url: $('#canvas_count_variations_id').attr("show_count_variations-url"),
	    success: function (data) {
	      
	      if (data['is_ok']) {
	    	  window.chartColors = {
	    			  red: 'rgb(255, 99, 132)',
	    			  orange: 'rgb(255, 159, 64)',
	    			  yellow: 'rgb(255, 205, 86)',
	    			  green: 'rgb(75, 192, 192)',
	    			  blue: 'rgb(54, 162, 235)',
	    			  purple: 'rgb(153, 102, 255)',
	    			  grey: 'rgb(231,233,237)'
		      };
	    	  
	    	  var color = Chart.helpers.color;
	    	  var barChartData = {
    	            labels: data['labels'],
    	            datasets: [{
    	                label: 'Variations <50%',
    	                borderColor: window.chartColors.red,
    	                backgroundColor: color(window.chartColors.red).alpha(0.5).rgbString(),
    	                borderWidth: 1,
    	                stack: 'Stack 0',
    	                data: data['data_less_50']
    	            }, {
    	                label: '50% < Variations < 90%',
    	                backgroundColor: color(window.chartColors.blue).alpha(0.5).rgbString(),
    	                borderColor: window.chartColors.blue,
    	                stack: 'Stack 0',
    	                data: data['data_50_var_90_50']
    	            }]
    	      };

	    	  // set the height
	    	  var height = data['sample_number'] * 6 + 30
	    	  $("#canvas_count_variations_id").attr("height", height) // 200 header
	    	  
	    	  /// load the chart
              var ctx = document.getElementById("canvas_count_variations_id").getContext("2d");
              window.myHorizontalBar = new Chart(ctx, {
                  type: 'horizontalBar',
                  data: barChartData,
                  options: {
                      // Elements options apply to all of the options unless overridden in a dataset
                      // In this case, we are setting the border of each horizontal bar to be 2px wide
                      elements: {
                          rectangle: {
                              borderWidth: 2,
                          }
                      },
                      responsive: true,
                      legend: {
                          position: 'right',
                      },
                      tooltips: {
                          mode: 'index',
                          intersect: false
                      },
                      title: {
                          display: true,
                          text: 'Variations <50% and 50%<Variations<90%'
                      },
                      scales: {
                          xAxes: [{
                              stacked: true,
                          }],
                          yAxes: [{
                              stacked: true
                          }]
                      }
                  }
              });
	      }
	      else{
	    	  $('#show_count_variations_id').append('<div class="alert alert-warning alert-dismissable"><strong>Fail</strong> ' + data['error_message'] + '</div>')
	      }
	    },
	    
	    // handle a non-successful response
	    error : function(xhr,errmsg,err) {
	    	$('#show_count_variations_id').append('<div class="alert alert-warning alert-dismissable"><strong>Fail</strong> to load the data from the server.</div>')
	        console.log(xhr.status + ": " + xhr.responseText); // provide a bit more info about the error to the console
	    }
    });
}


// show coverage graphic 
$(document).on("click", "a", function(){
	$('#modal-body-coverage').empty();
	
	$.ajax({
	  	beforeSend: function() {
	  		$('#loader_coverage_image').show();
	  	},
	  	complete: function(){
	  		$('#loader_coverage_image').hide();
	  	},
	
	  	data : { 
	  		project_sample_id : $(this).attr('project_sample_id'), // data sent with the get request 
	  		element : $(this).attr('sequence') }, // data sent with the get request 
	  		url: $('#modal-body-coverage').attr("show-coverage-modal-url"),
	  		success: function (data) {
	  			if (data['is_ok']) {
	  				$('h4.modal-title').text(data['text']);
	  				$('#modal-body-coverage').prepend(data['image'])
	  			}
	  		},
	      
	      // handle a non-successful response
	      error : function(xhr,errmsg,err) {
	          alert(errmsg);
	          console.log(xhr.status + ": " + xhr.responseText); // provide a bit more info about the error to the console
	      }
	});
});






