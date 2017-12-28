

$(function () {
  /* 1. OPEN THE FILE EXPLORER WINDOW */
  $(".js-upload-files").click(function () {
    $("#fileupload").click();
  });

  /* 2. INITIALIZE THE FILE UPLOAD COMPONENT */
  $("#fileupload").fileupload({
    dataType: 'json',
    done: function (e, data) {  /* 3. PROCESS THE RESPONSE FROM THE SERVER */
      if (data.result.is_valid) {
        $("#gallery tbody").prepend(
          "<tr><td><a href='" + data.result.url + "'>" + data.result.name + "</a></td>" +
             "<td>Uploaded</a></td>" +
          "</tr>"
        )
      }
      else{
    	  $("#gallery tbody").prepend(
    	          "<tr><td>" + data.result.name + "</a></td>" +
    	          	"<td>Failed: " + data.result.message + "</a></td>" +
    	          "</tr>"
          ) 
      }
    }
  });

});

/* Drag and drop files */
//$(function () {
//
//	  $(".js-upload-photos").click(function () {
//	    $("#fileupload").click();
//	  });
//
//	  $("#fileupload").fileupload({
//	    dataType: 'json',
//	    done: function (e, data) {
//	      if (data.result.is_valid) {
//	        $("#gallery tbody").prepend(
//	          "<tr><td><a href='" + data.result.url + "'>" + data.result.name + "</a></td></tr>"
//	        )
//	      }
//	    }
//	  });
//
//	});

/* Drag and drop files with progress bar */
//$(function () {
//
//	  $(".js-upload-photos").click(function () {
//	    $("#fileupload").click();
//	  });
//
//	  $("#fileupload").fileupload({
//	    dataType: 'json',
//	    sequentialUploads: true,
//
//	    start: function (e) {
//	      $("#modal-progress").modal("show");
//	    },
//
//	    stop: function (e) {
//	      $("#modal-progress").modal("hide");
//	    },
//
//	    progressall: function (e, data) {
//	      var progress = parseInt(data.loaded / data.total * 100, 10);
//	      var strProgress = progress + "%";
//	      $(".progress-bar").css({"width": strProgress});
//	      $(".progress-bar").text(strProgress);
//	    },
//
//	    done: function (e, data) {
//	      if (data.result.is_valid) {
//	        $("#gallery tbody").prepend(
//	          "<tr><td><a href='" + data.result.url + "'>" + data.result.name + "</a></td></tr>"
//	        )
//	      }
//	    }
//
//	  });
//
//	});