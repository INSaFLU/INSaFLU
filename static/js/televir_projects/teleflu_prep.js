
$(document).ready(function () {
  
  $('#teleflu_submit-button').click(function () {
      console.log('clicked');
      var url = $('#teleflu_submit-button').attr('href');
      var reload_url = $('#teleflu_submit-button').attr('reload_ref');
      var csrf = $('#teleflu_submit-button').attr('csrf');
      
      // get checked reference rows
      var checkedRows_refs = [];
      $('.teleflu_reference-checkbox:checked').each(function() {
        // collect ids of checked rows

        var ref_id= $(this).attr('ref_id');

        checkedRows_refs.push(ref_id);
      });
      // get checked samples rows
      var checkedRows_samples = [];
      $('.select_sample-checkbox:checked').each(function () {
        // collect ids of checked rows
        var sample_id = $(this).attr('sample_id');
        checkedRows_samples.push(sample_id);
      });

      console.log(checkedRows_samples);

      console.log(checkedRows_refs);
  
      // Process the checked rows
      // Add your processing logic here
      // send row ids to server using .ajax
      $.ajax({
        type: 'POST',
        url: url,
        data: {
          'csrfmiddlewaretoken': csrf,
          'ref_ids': checkedRows_refs,
          'sample_ids': checkedRows_samples,
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
              
          } else if (data['exists']) {
            $('#id_messages_remove').append('<div class="alert alert-dismissible alert-warning">' +
              'A Project with these references already exists.' +
              '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
              '</div>');

              // drop modal
              $('#myModal').modal('hide');
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
    var project_id = btn.attr("project_id");
    $.ajax({
        url: btn.attr("href"),
        type: 'GET',
        data: {
            search_add_project_reference: $('#teleflu_search-input').val(),
            teleflu_reference: "1",
            project_id: project_id,
            csrfmiddlewaretoken: csrf,
        },
        success: function (data) {

          $('#teleflu_reference_table_div').html(data.my_content);
          
          if (data["references_count"] > 0){
            
          var elements = document.getElementsByName("teleflu_select_ref");
          for(var i = 0, n = elements.length; i < n; i++){
            elements[i].addEventListener('click', toggle_reference_js_check_box, false);
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
          }
        }
    }
    )
  };
$('.load-teleflu-content').on('click', loadTelefluContent);
$('.teleflu_open').on('click', loadTelefluContent);

/// toggle checkbox
function toggle_reference_js_check_box(source) {
  var url = $('#teleflu_table_with_check_id').attr("set-check-box-values-url");
  var id = source.target.id;
  var elements = document.getElementsByName("teleflu_select_ref");
  for (var i = 0, n = elements.length; i < n; i++) {
    if (elements[i].id != id) {
      elements[i].checked = false;
    
    }
  }
};

/// toggle checkbox
function toggle_reference_check_box(source) {
  var url = $('#teleflu_table_with_check_id').attr("set-check-box-values-url");
  var id = source.target.id;
	$.ajax({
		url: url,
		data : { get_change_check_box_single : source.target.checked,
				 id : source.target.value }, // data sent with the post request
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