$('#request_map_selected').on("click", function(e){
    attr= $(this).attr('id');
    
    $('#id-label-map-selected').text('Map selected references?');
    
  });

  $('#request_screening').on("click", function(e){
    attr= $(this).attr('id');
    
    $('#id-label-screening').text('Deploy Screening?');
    
  });


  // deploy screening
  $('#id-screening-button').on('click', function() {

    var url= $('#id-modal-body-screening').attr("screening-single-value-url");
    var csrf_token = $('#id-modal-body-screening').attr('data-csrf'); 
    var sample_id = $('#headingsample').attr('sample-id');

    $.ajax({
      type: 'POST',
      url: url,
      data: {
        'csrfmiddlewaretoken': csrf_token,
        'sample_id': sample_id,
      },
      success: function(data) {
        if (data['is_ok'] == true && data['is_deployed'] == true) {
          $('#id_messages_remove').append('<div class="alert alert-dismissible alert-success">' +
            'Screening deployed. ' + data['message'] +
            '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
            '</div>');

        } else if (data['is_ok']) {
          $('#id_messages_remove').append('<div class="alert alert-dismissible alert-warning">' +
            'Screening not deployed. Check settings.' + data['message'] +
            '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
            '</div>');

            // drop modal
            $('#id_screening_modal').modal('hide');
          
        } else if (data['is_ok'] == false) {
          $('#id_messages_remove').append('<div class="alert alert-dismissible alert-warning">' +
            'Screening not deployed. Check Parameters.' +
            '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
            '</div>');

            // drop modal
            $('#id_screening_modal').modal('hide');
        }
      }
    })
  });

  $('.class-ref-checkbox').on('change', function() {
    var checkedRows = JSON.parse(localStorage.getItem('checkedRows')) || [];
    $('.class-ref-checkbox:checked').each(function () {
      
      var ref_id = $(this).attr('ref_id');
      if (!checkedRows.includes(ref_id)) {
        checkedRows.push(ref_id);
      }
    });
    // Store the array of checked rows in local storage
    localStorage.setItem('checkedRows', JSON.stringify(checkedRows));
  });

  function clearSelections() {
    // Assuming your checkboxes have a common class name 'table-checkbox'
    var checkedRows = []
    localStorage.setItem('checkedRows', JSON.stringify(checkedRows));
    $('.class-ref-checkbox').prop('checked', false);

  }


  $(document).ready(function() {
    var checkedRows = JSON.parse(localStorage.getItem('checkedRows')) || [];
    $('.class-ref-checkbox').each(function() {
        var ref_id = $(this).attr('ref_id');
        if (checkedRows.includes(ref_id)) {
            $(this).prop('checked', true);
        }
    });
    document.getElementById('clearSelections').addEventListener('click', clearSelections);
  });


  $('#id-map-selected-button').on('click', function() {

    var url= $('#id-modal-body-map-selected').attr("map-selected-single-value-url");
    var csrf_token = $('#id-modal-body-map-selected').attr('data-csrf'); 
    var sample_id = $('#headingsample').attr('sample-id');

    var checkedRows = JSON.parse(localStorage.getItem('checkedRows')) || [];

    $.ajax({
      type: 'POST',
      url: url,
      data: {
        'csrfmiddlewaretoken': csrf_token,
        'sample_id': sample_id,
        'reference_ids': checkedRows
      },
      success: function(data) {
        if (data['is_ok'] == true && data['is_deployed'] == true) {
          $('#id_messages_remove').append('<div class="alert alert-dismissible alert-success">' +
            'References mapping deployed. ' + data['message'] +
            '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
            '</div>');

        } else if (data['is_empty'] == true) {
          $('#id_messages_remove').append('<div class="alert alert-dismissible alert-warning">' +
            'No references were selected. ' + data['message'] +
            '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
            '</div>');

          // drop modal
          $('#id_map_selected_modal').modal('hide');
          
        } else if (data["is_already_mapped"] == true) {
          $('#id_messages_remove').append('<div class="alert alert-dismissible alert-warning">' +
            'Selected and added references already mapped.' +
            '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
            '</div>');
          
          // drop modal
          $('#id_map_selected_modal').modal('hide');
  
        } else if (data['is_ok']) {
          $('#id_messages_remove').append('<div class="alert alert-dismissible alert-warning">' +
            'References were not mapped. ' + data['message'] +
            '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
            '</div>');

            // drop modal
            $('#id_map_selected_modal').modal('hide');

        } else if (data['is_ok'] == false) {
          $('#id_messages_remove').append('<div class="alert alert-dismissible alert-warning">' +
            'The references were not mapped.' +
            '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
            '</div>');

            // drop modal
            $('#id_map_selected_modal').modal('hide');
        } 
        localStorage.removeItem('checkedRows');
        $('.class-ref-checkbox').prop('checked', false);
      }
    })
  });


  $('#remove-panel-submit-button').click(function (){
    var url = $('#id-modal-body-remove-panel').attr('remove-panel-single-value-url');
    var csrf_token = $('#id-modal-body-remove-panel').attr('data-csrf');
    var panel_id = $('#id-modal-body-remove-panel').attr('data-panel-id');
    var sample_id = $('#headingsample').attr('sample-id');
    var load_url = $('#panel-list-section').attr('data-url');
    $.ajax({
      type: 'POST',
      url: url,
      data: {
        'csrfmiddlewaretoken': csrf_token,
        'panel_id': panel_id,
        'sample_id': sample_id
      },
      success: function(data) {
  
        if (data['is_ok']) {
          $('#id_messages_remove').append('<div class="alert alert-dismissible alert-success">' +
            'Panel successfully removed' +
            '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
            '</div>');
            load_panels(sample_id, load_url, '.panel-list ul', true, false);
            // drop modal
            $('#removePanelModal').modal('hide');
        } else if (data['is_error']) {
          $('#id_messages_remove').append('<div class="alert alert-dismissible alert-warning">' +
            'The panel was not removed.' +
            '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
            '</div>');
            // drop modal
            $('#removePanelModal').modal('hide');
        }
      }
  });
  });


  $('#panel-submit-button').click(function () {
    var url = $(this).attr('href');
    var sample_id = $('#headingsample').attr('sample-id');
    var csrf_token = $('#headingsample').attr('data-csrf');
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
        'sample_id': sample_id,
        'panel_ids': checkedRows
      },
      success: function(data) {
        if (data['is_ok']) {
          $('#id_messages_remove').append('<div class="alert alert-dismissible alert-success">' +
            'Panels successfully added' +
            '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
            '</div>');
            load_panels(sample_id, load_url, '.panel-list ul', true, false);
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



  $('#submit-button').click(function() {
    var url = $('#submit-button').attr('href');
    var reload_url = $('#submit-button').attr('reload_ref');
    var sample_id = $('#headingsample').attr('sample-id');
    var csrf_token = $('#headingsample').attr('data-csrf');
    var checkedRows = [];
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
        'csrfmiddlewaretoken': csrf_token,
        'sample_id': sample_id,
        'reference_ids': checkedRows
      },
      success: function(data) {
        if (data['is_ok']) {
          $('#id_messages_remove').append('<div class="alert alert-dismissible alert-success">' +
            'References successfully added' +
            '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
            '</div>');
            // reload table
            update_ref_table();
            // clear selection 
            $('.reference-checkbox').prop('checked', false);
            
            
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


  $('#id-deploy-metagenomics-button').on('click', function(){
    url= $('#id-modal-body-deploy-metagenomics-sample').attr("deploy-metagenomics-single-value-url");
    sample_id = $('#id-modal-body-deploy-metagenomics-sample').attr('sample_id');
    csrf_token = $('#headingsample').attr('data-csrf');
  
    $.ajax({
        url: url, 
        type: 'POST',
        data: {
            sample_id : sample_id,
            csrfmiddlewaretoken: csrf_token,
        }, // data sent with the post request
        success: function(data) {

          
          if (data['is_ok'] === true && data['is_deployed'] === true){
              /// add message with informaton
              $('#id_messages_remove').append('<div class="alert alert-dismissible alert-success">' +
              'The sample \'' + $('#id-modal-body-deploy-metagenomics-sample').attr('ref_name') + '\' was successfully deployed.' +
              '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
              '</div>');
              
          }
          else if (data['is_ok'] === true && data['is_deployed'] === false){
              /// add message with informaton
              $('#id_messages_remove').append('<div class="alert alert-dismissible alert-warning">' +
              'The sample \'' + $('#id-modal-body-deploy-metagenomics-sample').attr('ref_name') + '\' was not deployed. Check if metagenomics settings are correctly set.' +
              '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
              '</div>');
          }
          else{
              /// add message with informaton
              $('#id_messages_remove').append('<div class="alert alert-dismissible alert-warning">' +
              'The sample \'' + $('#id-modal-body-deploy-metagenomics-sample').attr('ref_name') + '\' was not deployed. Check if metagenomics settings are correctly set.' +
              '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
              '</div>');
          }
        }
    })
  }
  );


  $('#request-map-panels').on('click', function() {
    var sample_name = $(this).attr('sample-name');
    var sample_id = $(this).attr('sample-id');
    $('#id-modal-body-map-panel').attr('sample_id', sample_id);
    $('#id-label-map-panel').text('Map sample \'' + sample_name + '\' to Added Reference Panels?');
  });

  $('#id-map-panel-button').on('click', function() {
    var sample_id = $('#id-modal-body-map-panel').attr('sample_id');
    var url = $('#id-modal-body-map-panel').attr('map-panel-single-value-url');
    var csrf = $('#id-modal-body-map-panel').attr('data-csrf');

    var data = {
      'sample_id': sample_id,
      'csrfmiddlewaretoken': csrf
    };
    $.ajax({
      type: 'POST',
      url: url,
      data: data,
      success: function(data) {
        if (data['is_ok'] === true && data['is_deployed'] === true) {
          $('#id_messages_remove').append('<div class="alert alert-success alert-dismissible" role="alert"><button type="button" ' +
            'class="close" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>'+
            'Panel mapping deployed successfully</div>');
        } else if (data['is_empty'] === true) { 
          $('#id_messages_remove').append('<div class="alert alert-warning alert-dismissible" role="alert"><button type="button" ' +
            'class="close" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>'+
            'No panels to map</div>');
        } else {
          $('#id_messages_remove').append('<div class="alert alert-danger alert-dismissible" role="alert"><button type="button" ' +
            'class="close" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>'+
            'Error deploying panel mapping</div>');
        }
      },
    }
    );
  });