
$(".remove_file").on('click', function(){

    var file_name = $(this).attr('file-name');
    var file_pk = $(this).attr('pk');
    var remove_single_value_url = $(this).attr('remove-single-value-url');
    $("#id-modal-body-remove-sample").html('<label id="id-label-remove" class="col-form-label">Remove file: '+file_name+'?</label>');
    $("#id-modal-body-remove-sample").attr('file_pk', file_pk);
    $("#id-remove-button").attr('remove-single-value-url', remove_single_value_url);
    $("#id_remove_modal").modal('show');
});

$("#id-remove-button").on('click', function(){
    var remove_single_value_url = $(this).attr('remove-single-value-url');
    var file_pk = $("#id-modal-body-remove-sample").attr('file_pk');

    $.ajax({
        url: remove_single_value_url,
        type: 'POST',
        data: {
            'file_id': file_pk,
            'csrfmiddlewaretoken': '{{ csrf_token }}'
        },
        success: function(data){
            if (data["is_ok"] == true){
                location.reload();
            } else {
                alert(data["message"]);
            }
        }
    });
});


function lockInput(elementId) {
    var element = $('#' + elementId);
    element.on('click keydown dragover drop', function(event) {
        // Prevent the default behavior for click, keydown, dragover, and drop events
        event.preventDefault();
    });
    // Add the 'locked-input' class to visually indicate the element is locked
    element.addClass('locked-input');
}

$(document).ready(function() {
    $('button[type="submit"]').prop('disabled', true);

    let url = $('#check-button').attr('url');
    
    $('#check-button').click(function() {
        var formData = new FormData($('#upload-form')[0]);
        $.ajax({
            url: url,
            type: 'POST',
            data: formData,
            processData: false,
            contentType: false,
            success: function(data) {
                $('#summary').html(data);
                if (data["is_ok"] == true) {
                    if (data["pass"] == true) {
                        $('button[type="submit"]').prop('disabled', false);

                        lockInput('fasta_file');
                        lockInput('metadata');
                    }
                    $('#summary').html(data["log"].join('<br>'));
                    $('#summary').addClass('active');
                    $('#summary').html(data["log"].join('<br>')).show();
                    // remove height of the summary div
                    $('#summary').css('height', 'auto');                  
                } else if (data["is_error"]) {
                    $('#id_messages_remove').append('<div class="alert alert-danger" role="alert">' + data["error_message"] + '</div>');
                } else {
                    $('button[type="submit"]').prop('disabled', true);
                }
            }
        });
    });
});
