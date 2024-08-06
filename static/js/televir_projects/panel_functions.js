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
                    icon.className = 'fa fa-search-plus';
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

                    var button = document.createElement("button");
                    button.className = "select-file-button btn btn-primary";
                    button.setAttribute("data-toggle", "modal");
                    button.setAttribute("data-target", "#panelSelectModal");
                    button.setAttribute
                    button.title = "Add from File";
                    
                    // Create the SVG
                    var svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
                    svg.setAttribute("class", "icon-white");
                    svg.setAttribute("version", "1.1");
                    svg.setAttribute("id", "Layer_1");
                    svg.setAttribute("xmlns", "http://www.w3.org/2000/svg");
                    svg.setAttribute("xmlns:xlink", "http://www.w3.org/1999/xlink");
                    svg.setAttribute("width", "15px");
                    svg.setAttribute("x", "0px");
                    svg.setAttribute("y", "0px");
                    svg.setAttribute("viewBox", "0 0 114.066 122.881");
                    svg.setAttribute("enable-background", "new 0 0 114.066 122.881");
                    svg.setAttribute("xml:space", "preserve");
                    
                    // Create the path
                    var path = document.createElementNS("http://www.w3.org/2000/svg", "path");
                    path.setAttribute("fill-rule", "evenodd");
                    path.setAttribute("clip-rule", "evenodd");
                    path.setAttribute("d", "M65.959,67.42h38.739c5.154,0,9.368,4.219,9.368,9.367v36.725 c0,5.154-4.221,9.369-9.368,9.369H65.959c-5.154,0-9.369-4.215-9.369-9.369V76.787C56.59,71.639,60.805,67.42,65.959,67.42 L65.959,67.42L65.959,67.42z M20.464,67.578c-1.495,0-2.74-1.352-2.74-2.988c0-1.672,1.209-2.989,2.74-2.989H43.88 c1.495,0,2.741,1.353,2.741,2.989c0,1.672-1.21,2.988-2.741,2.988H20.464L20.464,67.578L20.464,67.578z M87.795,18.186h9.822 c1.923,0,3.703,0.783,4.947,2.063c1.28,1.281,2.064,3.025,2.064,4.947v33.183h-6.051V25.196c0-0.285-0.107-0.533-0.285-0.711 c-0.177-0.178-0.426-0.285-0.711-0.285H87.76v34.18h-6.014V7.011c0-0.285-0.107-0.534-0.285-0.711 c-0.178-0.178-0.428-0.285-0.712-0.285H6.976c-0.285,0-0.535,0.106-0.712,0.285C6.085,6.478,5.979,6.726,5.979,7.011v83.348 c0,0.285,0.107,0.533,0.285,0.711s0.427,0.285,0.711,0.285h38.871v6.014H22.812v11.174c0,0.285,0.107,0.535,0.285,0.713 c0.177,0.176,0.427,0.285,0.711,0.285l22.038-0.002v6.014H23.844c-1.922,0-3.701-0.783-4.946-2.064 c-1.282-1.279-2.064-3.023-2.064-4.947l0-11.172H7.011c-1.922,0-3.701-0.785-4.946-2.064C0.783,94.023,0,92.279,0,90.357V7.011 C0,5.089,0.783,3.31,2.064,2.064C3.345,0.783,5.089,0,7.011,0h73.774c1.921,0,3.701,0.783,4.947,2.063 c1.28,1.282,2.063,3.025,2.063,4.947V18.186L87.795,18.186L87.795,18.186L87.795,18.186z M20.428,28.647 c-1.495,0-2.74-1.353-2.74-2.99c0-1.672,1.21-2.989,2.74-2.989l46.833,0c1.495,0,2.739,1.353,2.739,2.989 c0,1.672-1.208,2.99-2.739,2.99L20.428,28.647L20.428,28.647L20.428,28.647z M20.428,48.114c-1.495,0-2.74-1.353-2.74-2.989 c0-1.672,1.21-2.989,2.74-2.989l46.833,0c1.495,0,2.739,1.352,2.739,2.989c0,1.672-1.208,2.989-2.739,2.989L20.428,48.114 L20.428,48.114L20.428,48.114z M73.868,98.787c-2.007,0-3.634-1.627-3.634-3.635c0-2.006,1.627-3.633,3.634-3.633h7.823v-7.83 c0-2.006,1.628-3.633,3.634-3.633s3.634,1.627,3.634,3.633v7.83h7.829c2.007,0,3.634,1.627,3.634,3.633 c0,2.008-1.627,3.635-3.634,3.635h-7.823v7.83c0,2.006-1.627,3.633-3.634,3.633c-2.006,0-3.634-1.627-3.634-3.633v-7.83H73.868 L73.868,98.787L73.868,98.787z");
                    
                    // Append the path to the SVG
                    svg.appendChild(path);
                    
                    // Append the SVG to the button
                    button.appendChild(svg);
                    

                    div.appendChild(a);

                    div.appendChild(removeButton);
                    div.appendChild(addButton);
                    div.appendChild(button);
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

    $(".select-file-button").click(function(){
        var user_id = $("#new-panel-button").attr('user_id');
        var panel_id = $(this).closest('.panel-container').attr('data-panel-id');
        var url_user_files = $(".panel-list").attr('url-user-files');

        $.ajax({
            url: url_user_files,
            method: 'POST',
            data : {
                'user_id': user_id,
                'csrfmiddlewaretoken': document.querySelector('[name=csrfmiddlewaretoken]').value
            },

            success: function(data) {
                var html = '';
                for(var index in data.files) {
                    html += '<input type="radio" file-id="'+index+'" panel-id="'+panel_id+'" name="file" value="'+data.files[index]+'">';
                    html += '<label for="file'+index+'">'+data.files[index]+'</label><br>';
                }
                $("#file-select").html(html);
            }
        });
    });

}


var ready_document = function (user_id, reload_url) {

    $("#submit-file").click(function(){
        var user_id = $("#new-panel-button").attr('user_id');
        var selected_file_id= $('input[name="file"]:checked').attr('file-id');
        var panel_id = $('input[name="file"]:checked').attr('panel-id');
        var url_register_file_panel = $(this).attr('url-register-file-panel');

        console.log("selecting file");

        if (selected_file_id) {
            // Submit the selected file
            $.ajax({
                url: url_register_file_panel,
                type: 'post',
                data: {
                    'file_id': selected_file_id,
                    'panel_id': panel_id,
                    'csrfmiddlewaretoken': document.querySelector('[name=csrfmiddlewaretoken]').value
                },
                success: function(response) {
                    // Handle the success response
                    // For example, you can close the modal and show a success message
                    $('#panelSelectModal').modal('hide');
                    load_panel_refs(panel_id);
                    $('.modal-backdrop').remove();
                    reload_panels(user_id);
                    alert('File registered successfully');
                },
                error: function(response) {
                    // Handle the error response
                    // For example, you can show an error message
                    alert('Error registering file');
                }
            });
        } else {
            alert("Please select a file");
        }
    });

    $('.remove-panel-button').click(function () {

        var panel_id = $(this).attr('data-panel-id');
        // set panel_id to the remove button
        $('#remove-panel-button').attr('panel_id', panel_id);
        
    });

    $('#remove-panel-button').click(function () {
        var panel_id = $(this).attr('panel_id');
        var url = $(this).attr('url');
        var csrf = $(this).attr('csrf');
        
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
        
        if (checkedRows.length == 0) {
            $('#id_messages_remove').append('<div class="alert alert-dismissible alert-warning">' +
                'No references were selected.' +
                '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
                '</div>');
            return;
        }
  
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

            } else if (data['is_full']) {
                $('#id_messages_remove').append('<div class="alert alert-dismissible alert-warning">' +
                    'Some references were not added. The panel is full. The maximum number of references per panel is ' + data["max_references"] +
                    '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
                    '</div>');
                // drop modal
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
                $('body').removeClass('modal-open');
                $('.modal-backdrop').remove();
                // unblock the UI
                $.unblockUI();
                
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
}