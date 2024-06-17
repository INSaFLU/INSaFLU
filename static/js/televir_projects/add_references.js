
$(document).ready(function() {
    $('#submit-button').click(function () {
      var url = $('#submit-button').attr('href');
      var reload_url = $('#submit-button').attr('reload_ref');
      var checkedRows = [];
      var csrf = $('#submit-button').attr('csrf');
      var ref_index = $('#submit-button').attr('ref_index');
      $('.reference-checkbox:checked').each(function() {
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
              $('#reference_table_div').html(data.empty_content);
              // empty #search-input
              $('#search-input').val('');
              // clear selection 
              $('.reference-checkbox').prop('checked', false);
              
              
              // drop modal
              $('#myModal').modal('hide');
              console.log('reload_url', reload_url);

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
          $.unblockUI();
        }
      })

    });
});


var loadThisContent = function (event) {
    event.preventDefault();
    var btn = $(this);
    var csrf = btn.attr("csrf");
    $.ajax({
        url: btn.attr("href"),
        type: 'GET',
        data: {
            search_add_project_reference: $('#search-input').val(),
            csrfmiddlewaretoken: csrf,
        },
        success: function (data) {

          $('#reference_table_div').html(data.my_content);
          
          if (data["references_count"] > 0){

          var selectAllCheckbox = document.getElementById('Select-All-Checkbox-Modal');

          var referenceCheckboxes = document.getElementsByClassName('reference-checkbox');
     
          selectAllCheckbox.addEventListener('change', function() {
            for (var i = 0; i < referenceCheckboxes.length; i++) {
              referenceCheckboxes[i].checked = this.checked;
            }
          });
          }
        }
    }
    )
  };
  $('.load-content').on('click', loadThisContent);