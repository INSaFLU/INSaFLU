

$(function () {
  /* 1. OPEN THE FILE EXPLORER WINDOW */
  $(".js-upload-files").click(function () {
    $("#fileupload").click();
  });

  /* 2. INITIALIZE THE FILE UPLOAD COMPONENT */
  $("#fileupload").fileupload({
    dataType: 'json',
    sequentialUploads: true,  /* 1. SEND THE FILES ONE BY ONE */
    start: function (e) {  /* 2. WHEN THE UPLOADING PROCESS STARTS, SHOW THE MODAL */
        $("#modal-progress").modal("show");
    },
    stop: function (e) {  /* 3. WHEN THE UPLOADING PROCESS FINALIZE, HIDE THE MODAL */
    	window.setTimeout(function (){$('#modal-progress').modal("hide"); }, 500);
    },
    progressall: function (e, data) {  /* 4. UPDATE THE PROGRESS BAR */
        var progress = parseInt(data.loaded / data.total * 100, 10);
        var strProgress = progress + "%";
        $(".progress-bar").css({"width": strProgress});
        $(".progress-bar").text(strProgress);
    },
    done: function (e, data) {  /* 3. PROCESS THE RESPONSE FROM THE SERVER */
      if (data.result.is_valid) {
        $("#gallery tbody").prepend(
          "<tr><td><a href='" + data.result.url + "' download>" + data.result.name + "</a></td>" +
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
