var load_panels = function(sample_id, load_url, target, load = false, suggest= true) {
    // var url = '{% url "panel_list" %}';
    var panelList = document.querySelector(target);
    console.log('Loading panels for sample:', sample_id, 'from:', load_url, 'to:', panelList, 'load:', load, 'suggest:', suggest);
    $.ajax({
        url: load_url,
        method: 'GET',
        data: {
            'sample_id': sample_id
        },
        success: function(data) {
            // Clear the panel list
            
            while (panelList.firstChild) {
                panelList.removeChild(panelList.firstChild);
            }

            // Check if the panels array is empty
            if (data.panels.length === 0) {
                var li = document.createElement('li');
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

            // reload the connects            
            // clear messages 
            prep_loads();
            //$('#id_messages_remove').empty();
            }
        },
        error: function(error) {
            console.error('Error:', error);
        }
    });
  }

  var prep_loads= function() {
    $('.remove-panel-button').click(function () {
      var panel_id = $(this).attr('data-panel-id');
      $('#id-modal-body-remove-panel').attr('data-panel-id', panel_id);
      $('#id-label-remove-panel').text('Remove panel ' + panel_id + '?');
    });
  }

  var update_ref_table = function (event) {
    var sample_id = $('#headingsample').attr('sample-id');
    var reload_url = $('#submit-button').attr('reload_ref');
    var csrf_token = $('#headingsample').attr('data-csrf');

    $.ajax({
      url: reload_url,
      type: 'GET',
      data: {
          sample_id: sample_id,
          csrfmiddlewaretoken: csrf_token,
      },
      success: function (data) {

        $('#collapsehead2').html(data.my_content);
        $('#reference_table_div').html(data.empty_content);
        $('#add_references_title').html(
          "<strong>" + data.references_count + " Added References </strong>"
        );
        // empty #search-input
        $('#search-input').val('');
      }
    }
    );
  };

  var loadThisContent = function (event) {
    event.preventDefault();
    var btn = $(this);
    var url = btn.attr("href");
    var csrf_token = $('#headingsample').attr('data-csrf');

    $.ajax({
        url: btn.attr("href"),
        type: 'GET',
        data: {
            search_add_project_reference: $('#search-input').val(),
            csrfmiddlewaretoken: csrf_token,
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