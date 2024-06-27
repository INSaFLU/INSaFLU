var load_panels_main = function(load_url, target, load = false, suggest= true) {
    // var url = '{% url "panel_list" %}';
    var panelList = document.querySelector(target);
    $.ajax({
        url: load_url,
        method: 'GET',
        data: {},
        success: function(data) {
            // Clear the panel list
            while (panelList.firstChild) {
                panelList.removeChild(panelList.firstChild);
            }

            // Check if the panels array is empty
            if (data.panels.length === 0) {
                var li = document.createElement('li');
                li.className = 'no-panels';
                li.textContent = 'No panels available.';
                panelList.appendChild(li);

            } else {
                // Add the new panels to the list
                data.panels.forEach(function(panel) {
                    var li = document.createElement('li');
                    // div 

                    var div = document.createElement('div');
                    div.className = 'panel-container clearfix';
                    div.setAttribute('data-panel-id', panel.id);

                    var a = document.createElement('a');
                    a.className = 'panel-link';
                    a.href = '#';
                    a.setAttribute('data-panel-id', panel.id);
                    a.textContent = panel.name;

                    a.appendChild(document.createTextNode(' '));
                    var panel_icon = document.createElement('i');
                    panel_icon.className = 'fa ' + panel.icon;
                    a.appendChild(panel_icon);

                    
                    var addButton = document.createElement('button');
                    addButton.className = 'add-reference-button btn btn-primary';
                    addButton.setAttribute('data-panel-id', panel.id);
                    addButton.setAttribute('data-toggle', 'modal');
                    addButton.setAttribute('data-target', '#myModal');
                    var icon = document.createElement('i');
                    icon.className = 'fa fa-plus';
                    addButton.appendChild(icon);
                    
                    var removeButton = document.createElement('button');
                    removeButton.className = 'remove-panel-button btn btn-danger';
                    removeButton.setAttribute('data-panel-id', panel.id);
                    removeButton.setAttribute('data-toggle', 'modal');
                    removeButton.setAttribute('data-target', '#removePanelModal');
                    var icon = document.createElement('i');
                    icon.className = 'fa fa-trash';
                    removeButton.appendChild(icon);
                    
                    var refnumber_span = document.createElement('span');
                    refnumber_span.className = 'reference-note';
                    refnumber_span.textContent = panel.references_count;
                    refnumber_span.setAttribute('data-panel-id', panel.id);

                    var panel_checkbox = document.createElement('input');
                    panel_checkbox.type = 'checkbox';
                    panel_checkbox.className = 'panel-checkbox';
                    panel_checkbox.setAttribute('data-panel-id', panel.id);

                    if (load){
                      div.appendChild(removeButton);
                    }
                                        
                    div.appendChild(a);

                    if (suggest){
                      div.appendChild(panel_checkbox);
                    }
                    div.appendChild(refnumber_span);

                    li.appendChild(div);
                    panelList.appendChild(li);
                    
                });
            }
        },
        error: function(error) {
            console.error('Error:', error);
        }
    });
  }

  /// set wait screen
  $(".request-add-teleflu-sample").on("click", function () {
    var teleflu_id = $(this).attr('teleflu-id');
    $("#id-add-teleflu-sample-button").attr('teleflu-id', teleflu_id);
    
    var checkedRows_samples = [];
    $('.select_sample-checkbox:checked').each(function () {
        // collect ids of checked rows
        var sample_id = $(this).attr('sample_id');
        checkedRows_samples.push(sample_id);
    });
    // change text
    if (checkedRows_samples.length == 0) {
      $("#id-label-add-teleflu-sample").text("No samples selected.");
    } else {
      $("#id-label-add-teleflu-sample").text("Add selected samples ?");
    }
  });

  $("#id-add-teleflu-sample-button").on("click", function () {

      var teleflu_id = $(this).attr('teleflu-id');
      var url = $('#id-modal-body-add-teleflu-sample').attr('add-teleflu-single-value-url');
      var csrf = $('#teleflu_create-button').attr('csrf');

      // get checked samples rows
      var checkedRows_samples = [];
      $('.select_sample-checkbox:checked').each(function () {
          // collect ids of checked rows
          var sample_id = $(this).attr('sample_id');
          checkedRows_samples.push(sample_id);
      });
      $.ajax({
          url: url,
          type: 'POST',
          data: {
              'teleflu_id': teleflu_id,
              'sample_ids': checkedRows_samples,
              'csrfmiddlewaretoken': csrf
          },
          success: function (data) {
              if (data['not_added'] == true) {
                  alert('No samples added. They may already be in the project.');
              } else if (data['is_empty'] == true) {
                  alert('No samples selected.');
              } else {
                  alert('Samples added.');
                  location.reload();
              }
          }
      });

  });


  $('#panel-submit-button').click(function () {
    var url = $(this).attr('href');
    var project_id = $(this).attr('project-id');
    var csrf_token = $(this).attr('data-csrf');
    var checkedRows = [];
    var load_url = $('#panel-list-section').attr('data-url');
    $('.panel-checkbox:checked').each(function() {
      // collect ids of checked rows
      var panel_id = $(this).attr('data-panel-id');
      checkedRows.push(panel_id);
    });

    if (checkedRows.length === 0) {
      $('#id_messages_remove').append('<div class="alert alert-dismissible alert-warning">' +
        'No panels were selected.' +
        '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
        '</div>');
      return;
    }

    $.ajax({
      type: 'POST',
      url: url,
      data: {
        'csrfmiddlewaretoken': csrf_token,
        'project_id': project_id,
        'panel_ids': checkedRows
      },
      success: function(data) {
        if (data['is_ok']) {
          $('#id_messages_remove').append('<div class="alert alert-dismissible alert-success">' +
            'Panels successfully added' +
            '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
            '</div>');
            // clear selection 
            $('.panel-checkbox').prop('checked', false);
            // drop modal
            $('#addPanelModal').modal('hide');
        } else if (data['is_error']) {
          $('#id_messages_remove').append('<div class="alert alert-dismissible alert-warning">' +
            'The panels were not added.' +
            '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
            '</div>');
            // drop modal
            $('#addPanelModal').modal('hide');
        } else if (data['is_empty']) {
          $('#id_messages_remove').append('<div class="alert alert-dismissible alert-warning">' +
            'No panels were selected.' +
            '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
            '</div>');
            // drop modal
            $('#addPanelModal').modal('hide');
        }
      }
    });
  });


