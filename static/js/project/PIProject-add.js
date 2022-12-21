//
//	This test on the name change
//
$("#id_project_name").on("change paste keyup", function () {
    $.ajax({
      url: $("#id_form_table_project_reference").attr("data-validate-project-reference-url"),
     // data: form.serialize(),
      // dataType: 'json',
      data : { 
          project_name : $('#id_project_name').val(),
          csrfmiddlewaretoken: '{{ csrf_token }}' }, // data sent with the post request
      success: function (data) {
        $(document).find('#error_1_id_name').remove();
        if (data.is_taken) {
          var error_message = "<strong>" + data.error_message + "</strong>";
           $("#id_project_name_error").append(
              $('<span/>', { 										// creates a dynamic div element on the fly
              id: 'error_1_id_name',
              class: 'fields_error',
              html: error_message,
           }));
        }
      },
      
      // handle a non-successful response
      error : function(xhr,errmsg,err) {
          console.log(xhr.status + ": " + xhr.responseText); // provide a bit more info about the error to the console
      }
 });
});
