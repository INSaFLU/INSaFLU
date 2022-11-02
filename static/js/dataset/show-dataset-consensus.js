//
//	This test on the name change
//
//$("#phylocanvas").on("change paste keyup", function () {


// show coverage graphic 
$(document).on("click", "a", function(e){
	var attr = $(this).attr('id');
	
	if (typeof attr === 'undefined'){
		return;
	}
	
	/// remove samples
	if (attr === 'id_remove_consensus_modal'){
		
		var ref_name = $(this).attr('ref_name');
		var ref_project = $(this).attr('ref_project');
		// ID to remove line
		var tr_to_remove = 'row_' + $(this).attr('pk');
		
		$('#id-label-remove').text('Do you want to remove \'' + ref_name + '\'?');
		$('#id-modal-body-remove-sample').attr('pk', $(this).attr('pk'));
		$('#id-modal-body-remove-sample').attr('ref_name', ref_name);
		$('#id-modal-body-remove-sample').attr('ref_project', ref_project);
		$('#id-modal-body-remove-sample').attr('tr_to_remove', tr_to_remove);
	}
});


$('#collapseTwo').on('shown.bs.collapse', function () {
	draw_phylo_canvas();
});

//remove the tree from the screen
$('#collapseTwo').on('hidden.bs.collapse', function () {
	$('#phylocanvas').empty();
});

//draw phylocanvas
// set size of window "#phylocanvas" defined in 'static/css/flu-web-site.css'
function draw_phylo_canvas() {

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
	    	dataset_id : $('#phylocanvas').attr("dataset_id"),
			csrfmiddlewaretoken: '{{ csrf_token }}'
	    }, // data sent with the get request

	    headers : { "cache-control" : 'no-cache, must-revalidate, max-age=0' },
	    url: $('#phylocanvas').attr("show-phylo-canvas-url"),
	    success: function (data) {
	    	
	      if (data['is_ok']) {
	    	  (function (Phylocanvas) {

	    		  // create the tree
	    		  window.PhylocanvasTree = Phylocanvas.createTree('phylocanvas', {

	//	    		    history: {
	//	    		        parent: document.getElementById("map_phylocanvas"),
	//	    		        zIndex: -1000
	//	    		    },
		    			  contextMenu : [{
		    		            text: 'Normal Menu',
		    		            handler: 'triggerNormal',
		    		            internal: false,
		    		            leaf: false
		    		          }, {
		    		            text: 'Internal Menu',
		    		            handler: 'triggerInternal',
		    		            internal: true,
		    		            leaf: true
		    		          }, {
		    		            text: 'Save as PNG',
		    		            handler: 'exportCurrentTreeView',
		    		            internal: false,
		    		            leaf: false
		    		          }],
	    		  
		    		    history: false,
		    		    baseNodeSize: 2,
		    		    padding: 5,
		    		    font: "helvetica",
		    		    zoomFactor: 2,
		    		    labelPadding: 5,
		    		    showLabels: true,
		    		    alignLabels: true,
		    		    highlightSize: 1,
		    		    highlightWidth: -1,
		    		    fillCanvas: true // Fits hierarchical and rectangular trees.
		    		});
	    		  /// set the data
	    		  fetch(data['url_sample']).then(function(response) {
	    			  return response.json();
	    		  }).then(function(jsonData) {
	    			  phylTree(jsonData, data['tree'], data['root']);
	    			  window.PhylocanvasTree.load(data['tree']);
	    			  
	    			  window.PhylocanvasTree.setTreeType('circular'); // Choosing type of tree: takes radial, rectangular, circular, diagonal and hierarchy.
		    		  window.PhylocanvasTree.setNodeSize(document.getElementById("nodeSize").value); // Setting the node size accordingly to the value of the respective range element.
		    		  window.PhylocanvasTree.setTextSize(document.getElementById("textSize").value); // Setting labels' font size.

		    		  radioDetect(jsonData, Object.keys(Object.values(jsonData)[0]).length); // Show/hide individual metadata bars.
		    		  selectAlll(jsonData, Object.keys(Object.values(jsonData)[0]).length); // Show all metadata bars when respective button is clicked.
		    		  displayLabel(jsonData, Object.keys(Object.values(jsonData)[0]).length); // Show all metadata labels when respective button is clicked.
	    		  }).catch(err => {		/// if something goes wrong
	    			  $('#phylocanvas').empty();
	    			  $('#loader_phylocanvas').hide();
	    	    	  $('#phylocanvas').append('<div class="alert alert-warning alert-dismissable text-center"><strong>Fail</strong> to load the tree.</div>')
	    		  })

	    	 })(window.Phylocanvas);
	    	 $('#loader_phylocanvas').hide();
	    	 
	    	  //// set the files
			  $('#tree_nwk_id').empty()
			  $('#tree_nwk_id').append(data['tree_nwk_id'])
			  $('#tree_tree_id').empty()
			  $('#tree_tree_id').append(data['tree_tree_id'])
	      }
	      else{
	    	  $('#loader_phylocanvas').hide();
	    	  $('#phylocanvas').append('<div class="alert alert-warning alert-dismissable text-center"><strong>Fail</strong> to load the tree.</div>')
	      }
	    },
	    
	    // handle a non-successful response
	    error : function(xhr,errmsg,err) {
	    	$('#loader_phylocanvas').hide();
	    	$('#phylocanvas').append('<div class="alert alert-warning alert-dismissable text-center"><strong>Fail</strong> to load the tree.</div>')
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
	$('.smenubar').remove();
	$('#msa_viewer_nucleote_id').empty();
});

//draw nucleotide alignments
function draw_nucleotide_alignments() {
    
    $.ajax({
    	/// spin 
    	beforeSend: function() {
    		$('#msa_viewer_nucleote_id').empty();
    		$('.smenubar').remove();
    		$('#loader_msa_viewer_nucleote_id').show();
    	},
    	complete: function(){
    		$('#loader_msa_viewer_nucleote_id').hide();
    	},
    	
	    data : { 
	    	dataset_id : $('#msa_viewer_nucleote_id').attr("dataset_id"),
			csrfmiddlewaretoken: '{{ csrf_token }}'
	    }, // data sent with the get request
	    
	    url: $('#msa_viewer_nucleote_id').attr("show-msa-nucleotide-url"),
	    success: function (data) {
	      if (data['is_ok']) {
	    	  var rootDiv = document.getElementById("msa_viewer_nucleote_id");
	    	  var labelNameLength = data['max_length_label'] * 8 + 10;	// multiply by 13 because the labelFontsize is 13; https://github.com/wilzbach/msa in zoomer section
	    	  var m = msa({
	    			el: rootDiv,
	    			vis: {
	    				conserv: true,
	    				overviewbox: false,
	    			    labelId: false
	    			},
	    	        zoomer: {
	    	        	labelNameLength: labelNameLength,
	    	        },
	    	        menu: "small",
	    	        bootstrapMenu: true,
	    	  });
	    	  
	    	  msa.io.fasta.read(data['alignment_fasta_show_id'], function(err, seqs){
	    		  m.seqs.reset(seqs);
	    		  if (data.hasOwnProperty("gff3_show_id")){
		    		  d3.text(data['gff3_show_id'], function(datasetText) {
		    			  var features = gff.parseSeqs(datasetText, data['last_name_seq']);
		    			  m.seqs.addFeatures(features);
		    			  m.render();
		    		  });
	    		  }
	    		  else{
	    			  m.render();  
	    		  }
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
			$('#loader_msa_viewer_nucleote_id').hide();
	    	$('#msa_viewer_nucleote_id').append('<div class="alert alert-warning alert-dismissable"><strong>Fail</strong> to load the alignment.</div>')
	        console.log(xhr.status + ": " + xhr.responseText); // provide a bit more info about the error to the console
	    }
    });
}



$('#id-remove-button').on('click', function(){

	$.ajax({
        url: $('#id-modal-body-remove-sample').attr("remove-single-value-url"),
        data : { 
        	consensus_id : $('#id-modal-body-remove-sample').attr('pk'),
    		csrfmiddlewaretoken: '{{ csrf_token }}'
        }, // data sent with the post request
        		
        success: function (data) {
          if (data['is_ok']) {
        	  
        	  /// remove line
        	  document.getElementById($('#id-modal-body-remove-sample').attr('tr_to_remove')).remove();
        	  
        	  /// add message with informaton
        	  $('#id_messages_remove').append('<div class="alert alert-dismissible alert-success">' +
        		data['message'] +
				'<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
				'</div>');
          }
          else{
        	/// add message with informaton
        	  $('#id_messages_remove').append('<div class="alert alert-dismissible alert-warning">' +
        		data['message'] +
				'<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
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

