
$(document).ready(function() {
    $('#teleflu_submit-button').click(function () {
      var url = $('#teleflu_submit-button').attr('href');
      var reload_url = $('#teleflu_submit-button').attr('reload_ref');
      var checkedRows = [];
      var csrf = $('#teleflu_submit-button').attr('csrf');
      var ref_index = $('#teleflu_submit-button').attr('ref_index');
      $('.teleflu_reference-checkbox:checked').each(function() {
        // collect ids of checked rows

        var ref_id= $(this).attr('ref_id');

        checkedRows.push(ref_id);
      });
  
      // Process the checked rows
      // Add your processing logic here
      // send row ids to server using .ajax
      $.ajax({
        type: 'POST',
        url: url,
        data: {
          'csrfmiddlewaretoken': csrf,
          'ref_id': ref_index,
          'reference_ids': checkedRows,
        },
        success: function(data) {
          if (data['is_ok']) {
            $('#id_messages_remove').append('<div class="alert alert-dismissible alert-success">' +
              'References successfully added' +
              '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
              '</div>');
              // reload table
              $('#teleflu_reference_table_div').html(data.empty_content);
              // empty #search-input
              $('#teleflu_search-input').val('');
              // clear selection 
              $('.teleflu_reference-checkbox').prop('checked', false);
              
              
              // drop modal
              $('#myModal').modal('hide');

          } else if (data['is_error']) {
            $('#id_messages_remove').append('<div class="alert alert-dismissible alert-warning">' +
              'The references were not added.' +
              '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
              '</div>');

              // drop modal
              $('#myModal').modal('hide');
          } else if (data['is_empty']) {
            $('#id_messages_remove').append('<div class="alert alert-dismissible alert-warning">' +
              'No references were selected.' +
              '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
              '</div>');

              // drop modal
              
          }
          
        }
      })

    });
});


/// add function to toggle the checkbox in tables
/// set the Listener and get the checked in the server, and set the box check in the client
$(document).ready(function(){

  var elements = document.getElementsByName("teleflu_select_ref");

	for(var i = 0, n = elements.length; i < n; i++){
		elements[i].addEventListener('click', toggle_check_box, false);
  }
  
	// get if there's any checked in the server
	$.ajax({
		url: $('#teleflu_table_with_check_id').attr("set-check-box-values-url"),
		data : { 
			get_check_box_single : '1',
			csrfmiddlewaretoken: '{{ csrf_token }}' }, // data sent with the post request
		success: function (data) {
			for (key in data){
				if (key === 'is_ok'){ continue; }
				var remember = document.getElementById(key);
				if (remember != null){
					remember.checked = data[key];
				}
			}
		},
	});
});




var loadTelefluContent = function (event) {
    event.preventDefault();
    var btn = $(this);
    var csrf = btn.attr("csrf");
    $.ajax({
        url: btn.attr("href"),
        type: 'GET',
        data: {
            search_add_project_reference: $('#teleflu_search-input').val(),
            teleflu_reference: "1",
            csrfmiddlewaretoken: csrf,
        },
        success: function (data) {

          $('#teleflu_reference_table_div').html(data.my_content);
          
          if (data["references_count"] > 0){
            
          var elements = document.getElementsByName("teleflu_select_ref");
          console.log(elements);
          for(var i = 0, n = elements.length; i < n; i++){
            elements[i].addEventListener('click', toggle_check_box, false);
          }
          console.log("elements");

          // get if there's any checked in the server
          console.log($('#teleflu_table_with_check_id').attr("set-check-box-values-url"),);
          $.ajax({
            url: $('#teleflu_table_with_check_id').attr("set-check-box-values-url"),
            data : { 
              get_check_box_single : '1',
              csrfmiddlewaretoken: '{{ csrf_token }}' }, // data sent with the post request
            success: function (data) {
              for (key in data){
                if (key === 'is_ok'){ continue; }
                var remember = document.getElementById(key);
                if (remember != null){
                  remember.checked = data[key];
                }
              }
            },
          });
          }
        }
    }
    )
  };
$('.load-teleflu-content').on('click', loadTelefluContent);
  
/// toggle checkbox
function toggle_check_box(source) {
  console.log("toggle_check_box");
	$.ajax({
		url: $('#teleflu_table_with_check_id').attr("set-check-box-values-url"),
		data : { get_change_check_box_single : source.target.checked,
				 value : source.target.value }, // data sent with the post request
		success: function (data) { 
			for (key in data){
				if (key === 'is_ok'){ continue; }
				var remember = document.getElementById(key);
				if (remember != null){
					remember.checked = data[key];
				}
				else{
					remember.checked = false;
				}
			}
		},
	});
};