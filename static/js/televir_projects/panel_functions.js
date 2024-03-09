var reload_panels = function(userId) {
    // var url = '{% url "panel_list" %}';
    var reload_url = $("#create-panel-button").attr('reload-url');

    $.ajax({
        url: reload_url,
        method: 'GET',
        data: {
            'user_id': userId
        },
        success: function(data) {
            // Clear the panel list
            var panelList = document.querySelector('.panel-list ul');
            while (panelList.firstChild) {
                panelList.removeChild(panelList.firstChild);
            }
            var referenceList = document.querySelector('.reference-list ul');
            while (referenceList.firstChild) {
                referenceList.removeChild(referenceList.firstChild);
            }
            // Check if the panels array is empty
            if (data.panels.length === 0) {
                var li = document.createElement('li');
                li.textContent = 'No panels available.';
                panelList.appendChild(li);
                var li = document.createElement('li');
                li.textContent = 'No references available.';
                referenceList.appendChild(li);
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

                    
                    div.appendChild(a);
                    div.appendChild(removeButton);
                    div.appendChild(addButton);
                    div.appendChild(refnumber_span);

                    li.appendChild(div);
                    panelList.appendChild(li);
                    
                });

            load_panel_refs(data.panels[0].id);
            // reload the connects
            reload_connects();
            
            // clear messages 
            //$('#id_messages_remove').empty();
            }
        },
        error: function(error) {
            console.error('Error:', error);
        }
    });
}


var load_panel_refs = function (panelId) {
    var url = document.querySelector('.panel-list').getAttribute('data-url');
    $.ajax({
        url: url,
        method: 'GET',
        data: {
            'panel_id': panelId
        },
        success: function (data) {
            
            // change the reference-title
            var referenceTitle = document.querySelector('.reference-title');
            referenceTitle.textContent = 'References for panel: ' + data.panel_name;

            // Clear the references list
            var referenceList = document.querySelector('.reference-list ul');
            while (referenceList.firstChild) {
                referenceList.removeChild(referenceList.firstChild);
            }
        
            // Check if the references array is empty
            if (data.references.length === 0) {
                var li = document.createElement('li');
                li.textContent = 'No references available.';
                referenceList.appendChild(li);
            } else {
                // Add the new references to the list
                data.references.forEach(function(reference) {
                    var li = document.createElement('li');
                    li.className = 'reference-item clearfix';
        
                    var description = document.createElement('p');
                    description.className = 'reference-description'; 
                    description.textContent = reference.description;
                    li.appendChild(description);

                    var taxid = document.createElement('p');
                    taxid.className = 'reference-id'; 
                    taxid.textContent = 'Taxid: ' + reference.taxid;
                    li.appendChild(taxid);
        
                    var accid = document.createElement('p');
                    accid.className = 'reference-id'; 
                    accid.textContent = 'Accid: ' + reference.accid;
                    li.appendChild(accid);

                    var remove_button = document.createElement('button');
                    remove_button.className = 'remove-reference-button btn btn-danger';
                    remove_button.id = 'remove-reference-button';
                    remove_button.setAttribute('ref_id', reference.id);
                    remove_button.setAttribute('panel_id', panelId);
                    remove_button.setAttribute('data-toggle', 'modal');
                    remove_button.setAttribute('data-target', '#removeReferenceModal');
                    var icon = document.createElement('i');
                    icon.className = 'fa fa-trash';
                    remove_button.appendChild(icon);
                    li.appendChild(remove_button);
        
                    referenceList.appendChild(li);
                });
            }
            // reload the connects
            $('.remove-reference-button').click(function () {
                var ref_id = $(this).attr('ref_id');
                var panel_id = $(this).attr('panel_id');

                $('#removeReferenceModal').attr('ref_id', ref_id);
                $('#removeReferenceModal').attr('panel_id', panel_id);

            
            });

        },
        error: function(error) {
            console.error('Error:', error);
        }
    });
}    

var loadThisContent = function (event) {
    event.preventDefault();
    var btn = $(this);
    var csrf = btn.attr("csrf");
    $.ajax({
        url: btn.attr("href"),
        type: 'GET',
        data: {
            search_add_project_reference: $('#search-input').val(),
            max_references: 10,
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


var reload_connects = function () {


    $('.panel-container').on('click', function () {
        var panelId = $(this).data('panel-id');
        load_panel_refs(panelId);
        
    });
    
    $('.remove-panel-button').click(function (e) {
        // stop event propagation
        
        var panel_id = $(this).attr('data-panel-id');
        // set panel_id to the remove button
        $('#remove-panel-button').attr('panel_id', panel_id);
        
    });

    $('.add-reference-button').click(function () {
        var panel_id = $(this).attr('data-panel-id');
        $('#submit-button').attr('ref_index', panel_id);
    });

}


var ready_document = function (user_id, reload_url) {

    $('.remove-panel-button').click(function () {

        var panel_id = $(this).attr('data-panel-id');
        // set panel_id to the remove button
        $('#remove-panel-button').attr('panel_id', panel_id);
        
    });

    $('#remove-panel-button').click(function () {
        var panel_id = $(this).attr('panel_id');
        var url = $(this).attr('url');
        var csrf = $(this).attr('csrf');
        var reload_url = $("create-panel-button").attr('reload-url');
        
        $.ajax({
            type: 'POST',
            url: url,
            data: {
                'csrfmiddlewaretoken': csrf,
                'panel_id': panel_id,
            },
            success: function(data) {
                if (data['is_ok']) {
                    $('#id_messages_remove').append('<div class="alert alert-dismissible alert-success">' +
                        'Panel successfully removed' +
                        '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
                        '</div>');
                    // reload table
                    $('#reference_table_div').html(data.empty_content);
                    // empty #search-input
                    $('#search-input').val('');
                    // clear selection 
                    $('.reference-checkbox').prop('checked', false);
                    // drop modal
                    load_panel_refs(panel_id);
                    // close modal 
                    $('#removePanelModal').modal('hide');
                    $('body').removeClass('modal-open');
                    $('.modal-backdrop').remove();
                    // unblock the UI
                    reload_panels(user_id);
                    $.unblockUI();
                } else if (data['is_error']) {
                    $('#id_messages_remove').append('<div class="alert alert-dismissible alert-warning">' +
                        'The panel was not removed.' +
                        '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
                        '</div>');
                }
            }
        })
    })


    $('#submit-button').click(function () {
      var url = $('#submit-button').attr('href');
      //var reload_url = $('#submit-button').attr('reload_ref');
      
      var checkedRows = [];
      var csrf = $('#submit-button').attr('csrf');
      var panel_index = $('#submit-button').attr('ref_index');
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
          'ref_id': panel_index,
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
              load_panel_refs(panel_index);
              reload_panels(user_id);
              // close modal 
              $('#myModal').modal('hide');
              $('body').removeClass('modal-open');
              $('.modal-backdrop').remove();
              // unblock the UI
              $.unblockUI();


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
}