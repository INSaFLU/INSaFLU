
$(document).on("click", "a", function (e) {

    if ($(this).attr("id") === "id_set_control") {

        var attr = $(this).attr('id');
        var ref_name = $(this).attr('ref_name');
        var sample_pk = $(this).attr('pk');
        var csrf_token = $("#submit-button").attr('csrf');
        
        // For some browsers, `attr` is undefined; for others `attr` is false.  Check for both.

        // add show to class 
        //check if class icon is filled circle or empty circle
        var icon = $(this).find('i');
        var icon_class = icon.attr('class');

        
        if (icon_class === 'fa fa-circle-o'){
            $('#id-label-set-control').text('Set \'' + ref_name + '\' as control?');
        }
        else {
            $('#id-label-set-control').text('Unset \'' + ref_name + '\' as control?');
            $('#id-set-control-button').text('Unset');
        }
        
        $('#id-modal-body-set-control').attr('pk', sample_pk);
        $('#id-modal-body-set-control').attr('ref_name', ref_name);
        $('#id-modal-body-set-control').attr('token', csrf_token);
        
    }
    else if ($(this).attr("id") === 'id_add_set_control_message'){
        $('#id_messages_set_control').append('<div class="alert alert-dismissible alert-warning">' +
                'No jobs to deploy.' +
                '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
                '</div>');
    }
    else if ($(this).attr("id") === 'deploy_metagenomics_modal'){
        var ref_name = $(this).attr('ref_name');
        var sample_pk = $(this).attr('pk');

        $('#id-label-deploy-metagenomics').text('Deploy metagenomics combined analysis for sample \'' + ref_name + '\'?');
        $('#id-modal-body-deploy-metagenomics-sample').attr('sample_id', sample_pk);
        $('#id-modal-body-deploy-metagenomics-sample').attr('ref_name', ref_name);
    }
    
});

$('#id-deploy-metagenomics-button').on('click', function () {
    url = $('#id-modal-body-deploy-metagenomics-sample').attr("deploy-metagenomics-single-value-url");
    sample_id = $('#id-modal-body-deploy-metagenomics-sample').attr('sample_id');
    csrf_token = $('#teleflu_create-button').attr("csrf");

    $.ajax({
        url: url,
        type: 'POST',
        data: {
            sample_id: sample_id,
            csrfmiddlewaretoken: csrf_token,
        }, // data sent with the post request
        success: function (data) {
            if (data['is_ok'] === true && data['is_deployed'] === true) {
                /// add message with informaton
                $('#id_messages_remove').append('<div class="alert alert-dismissible alert-success">' +
                    'The sample \'' + $('#id-modal-body-deploy-metagenomics-sample').attr('ref_name') + '\' was successfully deployed.' +
                    '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
                    '</div>');
                location.reload()
            }
            else if (data['is_ok'] === true && data['is_deployed'] === false) {
                /// add message with informaton
                $('#id_messages_remove').append('<div class="alert alert-dismissible alert-warning">' +
                    'The sample \'' + $('#id-modal-body-deploy-metagenomics-sample').attr('ref_name') + '\' was not deployed. Check if metagenomics settings are correctly set.' +
                    '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
                    '</div>');
            }
            else {
                /// add message with informaton
                $('#id_messages_remove').append('<div class="alert alert-dismissible alert-warning">' +
                    'The sample \'' + $('#id-modal-body-deploy-metagenomics-sample').attr('ref_name') + '\' was not deployed. Check if metagenomics settings are correctly set.' +
                    '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
                    '</div>');
            }
        }
    })
});

$('.kill-button').on('click', function () {
    event.preventDefault();

    csrf_token = $('#teleflu_create-button').attr("csrf");
    url = $(this).attr("deploy-url");

    $.ajax({
        url: $(this).attr("deploy-url"),
        type: "POST",
        data: {
            'csrfmiddlewaretoken': csrf_token,
            'sample_id': $(this).attr("sample_id"),
        },
        data_type: 'json',
        success: function (data) {
            if (data["is_ok"] == true && data["is_deployed"] == false) {
                alert("No runs to deploy, check project settings.");
            }
            else if (data["is_ok"] == true && data["is_deployed"] == true) {
                alert("Runs deployed");
            }

        }

    });
});

$('#sort_sample_btn').on('click', function () {

    event.preventDefault();
        
    url = $(this).attr("sort-url");
    csrf_token = $('#teleflu_create-button').attr("csrf");

    $.ajax({
        url: $(this).attr("sort-url"),
        type: "POST",
        data: {
            'csrfmiddlewaretoken': csrf_token,
            'sample_id': $(this).attr("sample_id"),
        },
        data_type: 'json',
        success: function (data) {
            if (data["is_ok"] == true && data["is_deployed"] == false) {
                alert("No Results to Sort");
            }
            else if (data["is_ok"] == true && data["is_deployed"] == true) {
                alert("Sort deployed");
            }
        }

    });

});


$(document).on("click", "a", function (e) {

    /// set wait screen
    //var id_ = $(this).attr('id');
    var href = $(this).attr('href');
    var onclick = $(this).attr('onclick');
    var id_ = $(this).attr('id');
    // check if href is defined

    if (typeof href !== 'undefined' && onclick !== 'undefined' && id_ !== 'undefined') {
        
        if (href !== '#id_set_default_modal' && id_ !== "deploypi_sample_btn" && id_ !== "sort_sample_btn" && onclick !== 'return false;' && id_ !== 'sidenavToggler' &&
            !href.startsWith('/media') && !href.startsWith('http')) {
            wait_screen();
        };

    };

});


$("#deploypi_mapping_btn").click(function (e) {

    // get checked samples rows
    var checkedRows_samples = [];
    $('.select_sample-checkbox:checked').each(function () {
        // collect ids of checked rows
        var sample_id = $(this).attr('sample_id');
        checkedRows_samples.push(sample_id);
    });

    $.ajax({
        url: $('#deploypi_mapping_btn').attr("deploy-url"),
        type: "POST",
        data: {
            'csrfmiddlewaretoken': $('#teleflu_create-button').attr("csrf"),
            'id': $(this).attr('id'),
            'user_id': $('#deploypi_btn').attr('user-id'),
            'project_id': $('#teleflu_create-button').attr("ref_index"),
            'sample_ids': checkedRows_samples,
        },
        data_type: 'json',
        success: function (data) {
            if (data["is_ok"] == true && data["is_deployed"] == false) {
                alert("No runs to deploy, check project settings.");
            }
            else if (data["is_ok"] == true && data["is_deployed"] == true) {
                alert("Runs deployed");
            }
            $.unblockUI();
        }

    });
});

$("#deploypi_added_mapping_btn").click(function (e) {
    var user_id = $('#deploypi_btn').attr('user-id');
    var project_id = $('#deploypi_btn').attr('project-id');
    csrf_token = $('#teleflu_create-button').attr("csrf");

    // get checked samples rows
    var checkedRows_samples = [];
    $('.select_sample-checkbox:checked').each(function () {
        // collect ids of checked rows
        var sample_id = $(this).attr('sample_id');
        checkedRows_samples.push(sample_id);
    });
    
    $.ajax({
        url: $('#deploypi_added_mapping_btn').attr("deploy-url"),
        type: "POST",
        data: {
            'csrfmiddlewaretoken': csrf_token,
            'id': $(this).attr('id'),
            'user_id': user_id,
            'project_id': project_id,
            'sample_ids': checkedRows_samples,
        },
        data_type: 'json',
        success: function (data) {
            if (data["is_ok"] == true && data["is_deployed"] == false) {
                alert(data["message"]);
            }
            else if (data["is_ok"] == true && data["is_deployed"] == true) {
                var how_many = data["samples_deployed"];
                alert(data["message"]);
            }
            $.unblockUI();
        }

    });
});


$("#deploypi_panels_btn").click(function (e) {
    var user_id = $('#deploypi_btn').attr('user-id');
    var project_id = $('#deploypi_btn').attr('project-id');
    csrf_token = $('#teleflu_create-button').attr("csrf");

    // get checked samples rows
    var checkedRows_samples = [];
    $('.select_sample-checkbox:checked').each(function () {
        // collect ids of checked rows
        var sample_id = $(this).attr('sample_id');
        checkedRows_samples.push(sample_id);
    });
    
    $.ajax({
        url: $('#deploypi_panels_btn').attr("deploy-url"),
        type: "POST",
        data: {
            'csrfmiddlewaretoken': csrf_token,
            'id': $(this).attr('id'),
            'user_id': user_id,
            'project_id': project_id,
            'sample_ids': checkedRows_samples,
        },
        data_type: 'json',
        success: function (data) {
            if (data["is_ok"] == true && data["is_deployed"] == false) {
                alert(data["message"]);
            }
            else if (data["is_ok"] == true && data["is_deployed"] == true) {
                var how_many = data["samples_deployed"];
                alert(data["message"]);
            }
            $.unblockUI();
        }

    });
});


$("#deploypi_btn").click(function (e) {

    var user_id = $('#deploypi_btn').attr('user-id');
    var project_id = $('#deploypi_btn').attr('project-id');
    csrf_token = $('#teleflu_create-button').attr("csrf");

    // get checked samples rows
    var checkedRows_samples = [];
    $('.select_sample-checkbox:checked').each(function () {
        // collect ids of checked rows
        var sample_id = $(this).attr('sample_id');
        checkedRows_samples.push(sample_id);
    });

    $.ajax({
        url: $('#deploypi_btn').attr("deploy-url"),
        type: "POST",
        data: {
            'csrfmiddlewaretoken': csrf_token,
            'id': $(this).attr('id'),
            'user_id': user_id,
            'project_id': project_id,
            'sample_ids': checkedRows_samples,
        },
        data_type: 'json',
        success: function (data) {
            if (data["is_ok"] == true && data["is_deployed"] == false) {
                alert("No runs to deploy, check project settings.");
            }
            else if (data["is_ok"] == true && data["is_deployed"] == true) {
                alert("Runs deployed");
            }
            $.unblockUI();
        }

    });
});

$("#sortpi_btn").click(function () {

    var user_id = $('#deploypi_btn').attr('user-id');
    var project_id = $('#deploypi_btn').attr('project-id');
    csrf_token = $('#teleflu_create-button').attr("csrf");


    $.ajax({
        url: $('#sortpi_btn').attr("deploy-url"),
        type: "POST",
        data: {
            'csrfmiddlewaretoken': csrf_token,
            'id': $(this).attr('id'),
            'user_id': user_id,
            'project_id': project_id,
        },
        data_type: 'json',
        success: function (data) {
            if (data["is_ok"] == true && data["is_deployed"] == false) {
                alert("No sort to deploy");
            }
            else if (data["is_ok"] == true && data["is_deployed"] == true) {
                alert("Sort deployed");
            }
        }

    });
});

$('#confirm-delete-teleflu-project-button').click(function() {
    const projectId = $(this).attr('project-id'); // Assuming 'project-id' is stored as a data attribute
    const deleteUrl = $(this).attr('delete-url'); // Get the delete URL from the button's attribute

    // Make the AJAX call
    $.ajax({
        url: deleteUrl,
        type: 'POST', // Assuming the server expects a POST request for deletion
        data: {
            'project_id': projectId, // Send the project ID as part of the request
            'csrfmiddlewaretoken': $('input[name=csrfmiddlewaretoken]').val() // Include CSRF token for Django
        },
        success: function(data) {
            if (data['is_ok']) {
                alert('Project deleted successfully');
                teleflu_projects_load();
            }
        },
        error: function(xhr, errmsg, err) {
            // Handle error (e.g., show an error message)
            alert('Error deleting project: ' + errmsg);
        }
    });
  });