/// set sample control

$(document).on("click", "a", function (e) {
    if ($(this).attr("id") === "id_set_control_modal") {
        
        var attr = $(this).attr('id');
        var ref_name = $(this).attr('ref_name');
        var sample_pk = $(this).attr('pk');
        
        // For some browsers, `attr` is undefined; for others `attr` is false.  Check for both.
        if (attr === 'id_set_control_modal'){
            $('#id-label-set-control').text('Set \'' + ref_name + '\' as control?');
            $('#id-modal-body-set-control').attr('pk', sample_pk);
            $('#id-modal-body-set-control').attr('ref_name', ref_name);
        }
        else if (attr === 'id_add_set_control_message'){
            $('#id_messages_set_control').append('<div class="alert alert-dismissible alert-warning">' +
                    'No jobs to deploy.' +
                    '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
                    '</div>');
        }
    }
});


$('#id-set-control-button').on('click', function(){

    url = $('#id-modal-body-set-control').attr("set-control-url");
    sample_id = $('#id-modal-body-set-control').attr('pk');

    $.ajax({
        url: url,
        type: 'POST',
        data : {
            sample_id : sample_id,
            csrfmiddlewaretoken: token,
        }, // data sent with the post request
                
        success: function (data) {
          if (data['is_ok']) {
              
              /// add message with informaton
              $('#id_messages_set_control').append('<div class="alert alert-dismissible alert-success">' +
                'The sample \'' + $('#id-modal-body-set-control').attr('ref_name') + '\' was successfully set as control.' +
                '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
                '</div>');
            }
            else{
                /// add message with informaton
                $('#id_messages_set_control').append('<div class="alert alert-dismissible alert-warning">' +
                'The sample \'' + $('#id-modal-body-set-control').attr('ref_name') + '\' was not set as control.' +
                '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
                '</div>');
            }
        },

        error: function (data) {
            console.log('error');
            console.log(data);
        }
    });
});
