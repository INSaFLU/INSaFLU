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

  function createTelefluProjectButtons(project, panelContainer) {
    // Assuming 'project' is an object that includes 'id' and 'samples'

    // Create the 'View Details' button
    var viewDetailsButton = document.createElement('a');
    viewDetailsButton.href = `/pathogen_identification/teleflu_project/${project.id}`; // Adjust the URL pattern as needed
    viewDetailsButton.className = 'teleflu-details btn btn-primary';
    viewDetailsButton.textContent = 'View Details';

    // Create the 'Add Samples' button
    var addSamplesButton = document.createElement('a');
    addSamplesButton.setAttribute('teleflu-id', project.id);
    addSamplesButton.href = "#";
    addSamplesButton.className = 'request-add-teleflu-sample btn btn-primary';
    addSamplesButton.setAttribute('data-toggle', 'modal');
    addSamplesButton.setAttribute('data-target', '#add_teleflu_sample_modal');

    var icon = document.createElement('i');
    icon.className = 'fa fa-plus';
    icon.setAttribute('aria-hidden', 'true');
    icon.style.marginRight = '5px';
    addSamplesButton.appendChild(icon);

    var samplesText = document.createTextNode(` Samples: ${project.samples}`);
    addSamplesButton.appendChild(samplesText);

    // Assuming 'panelContainer' is the container to which you want to append the buttons
    panelContainer.append(viewDetailsButton);
    panelContainer.append(addSamplesButton);

    // Append 'panelContainer' to the desired parent element in your document
    // For example, document.getElementById('someParentElementId').appendChild(panelContainer);
  }


  function addInsafluProjectStatusToPanel(project, panelContainer) {
    // Create the 'teleflu-results' div
    var resultsDiv = document.createElement('div');
    resultsDiv.className = 'teleflu-results';

    // Check the project's 'insaflu_project' status and create the corresponding elements
    if (project.insaflu_project === "Finished") {
        var span = document.createElement('span');
        span.className = 'insaflu-exists';
        var a = document.createElement('a');
        var i = document.createElement('i');
        i.className = 'fa fa-check';
        i.setAttribute('aria-hidden', 'true');
        a.append(i);
        a.append(' INSaFLU Project');
        span.append(a);
        resultsDiv.append(span);
    } else if (project.insaflu_project === "Processing") {
        var span = document.createElement('span');
        span.className = 'insaflu-processing';
        var a = document.createElement('a');
        var i = document.createElement('i');
        i.className = 'fa fa-spinner fa-spin';
        i.setAttribute('aria-hidden', 'true');
        a.append(i);
        a.append(' INSaFLU Project');
        span.append(a);
        resultsDiv.append(span);
    }

    // Append the 'teleflu-results' div to the 'panelContainer'
    panelContainer.append(resultsDiv);
}

var teleflu_projects_load = function() {
  var url = $('#teleflu-projects-info').attr('teleflu-projects-url');

  $.ajax({
    url: url, // The URL to your endpoint that returns teleflu_projects data
    type: 'GET', // or 'POST', depending on your server setup
    dataType: 'json', // Expecting JSON data in response
    data: {
      'project_id': $('#panel-submit-button').attr('project-id'),
    },
    success: function (data) {

      if (data["is_ok"] === true && data["is_empty"] === false) {
          // Create the container div
          var containerDiv = $('<div>', { class: 'teleflu-table-container' });
          
          // Add title and description
          containerDiv.append('<h4 class="table-title">Reference Analysis</h4>');
          containerDiv.append('<p class="table-description"><i class="fa fa-crosshairs fa-lg" aria-hidden="true">‌</i>Focus on selected targets for validation. Coordinate processing and mapping workflows against selected references across multiple samples.</p>');
          
          // Create the list
          var projectsList = $('<ul>', { id: 'teleflu-list' });
          $.each(data.teleflu_projects, function(i, project) {
              var listItem = $('<li>');

              var panelContainer = $('<div>', { class: 'teleflu-panel-container' });
              var teleflu_manage = $('<div>', { class: 'teleflu-manage clearfix' });
            
              // Create the 'Delete Project' button
              var deleteProjectButton = document.createElement('button');
              deleteProjectButton.className = 'btn delete-project-button';
              deleteProjectButton.setAttribute('data-toggle', 'modal');
              deleteProjectButton.setAttribute('data-target', '#delete_teleflu_project_modal');
              // Optionally, set a custom attribute to hold the project ID for deletion
              deleteProjectButton.setAttribute('project-id', project.id);
              deleteProjectButton.setAttribute('data-project-name', project.ref_accid);
              teleflu_manage.append(deleteProjectButton); // Append the delete button
              
              // Add project details to panelContainer...
              teleflu_manage.append('<div class="reference-item"><p class="reference-description">' + project.ref_description + '</p><p class="reference-id">Taxid: ' + project.ref_taxid + '</p><p class="reference-id">Accid: ' + project.ref_accid + '</p></div>');
              
              // Add buttons to panelContainer
              createTelefluProjectButtons(project, teleflu_manage);            
              panelContainer.append(teleflu_manage);
            
              addInsafluProjectStatusToPanel(project, panelContainer);

              // Add more project details as needed...
              listItem.append(panelContainer);
              projectsList.append(listItem);
          });
          
          // Append the list to the container
          containerDiv.append(projectsList);
          // clear previous content
          $('#teleflu-projects-info').empty();
          // Append the container to the teleflu-projects-info div or any other target element
          $('#teleflu-projects-info').append(containerDiv);
          
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
        
          
          $(".delete-project-button").on("click", function () {
            var project_id = $(this).attr('project-id');
            var project_name = $(this).attr('data-project-name');

            // Update the modal's body to include the project name
            const modalBody = document.querySelector('#delete_teleflu_project_modal .modal-body');
            modalBody.innerHTML = `Are you sure you want to delete the project: <strong>${project_name}</strong>?`;

            // Set the 'project-id' attribute on the 'confirm-delete-teleflu-project-button'
            const confirmDeleteButton = document.querySelector('#confirm-delete-teleflu-project-button');
            confirmDeleteButton.setAttribute('project-id',project_id);
          });

      } else if (data["is_empty"] === true) {
          // Display a message when no teleflu projects are available
          $('#teleflu-projects-info').empty();
      }
    },
    error: function(xhr, status, error) {
        console.error("Error fetching teleflu projects: ", error);
    }
});
}
