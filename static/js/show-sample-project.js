//
//	This test on the name change
//
//$("#phylocanvas").on("change paste keyup", function () {
$('#collapseTwo').on('shown.bs.collapse', function () {
      $.ajax({
    	beforeSend: function() {
    		$('#loader_phylocanvas').show();
    	},
    	complete: function(){
    		$('#loader_phylocanvas').hide();
    	},

    	data : { project_id : $('#phylocanvas').attr("project_id") }, // data sent with the get request 
        url: $('#phylocanvas').attr("show-phylo-canvas-url"),
        success: function (data) {
          if (data['is_ok']) {
        	  (function (Phylocanvas) {
        	      var tree = Phylocanvas.createTree('phylocanvas', { history: false, });
        	      tree.load(data['tree']);
        	      tree.setTreeType('rectangular');
        	  })(window.Phylocanvas);
          }
        },
        
        // handle a non-successful response
        error : function(xhr,errmsg,err) {
            alert(errmsg);
            console.log(xhr.status + ": " + xhr.responseText); // provide a bit more info about the error to the console
        }
   });
});

// remove the tree from the screen
$('#collapseTwo').on('hidden.bs.collapse', function () {
	$('#phylocanvas').empty();
});

