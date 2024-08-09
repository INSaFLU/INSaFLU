

///
$(document).on("click", "a", function (e) {
  var attr = $(this).attr('id');
  var ref_name = $(this).attr('ref_name');
  var tr_to_remove = e.target.parentNode.parentNode.parentNode.id;

  // For some browsers, `attr` is undefined; for others `attr` is false.  Check for both.
  if (attr === 'id_remove_sample_modal') {
    $('#id-label-remove').text('Do you want to remove \'' + ref_name + '\'?');
    $('#id-modal-body-remove-sample').attr('pk', $(this).attr('pk'));
    $('#id-modal-body-remove-sample').attr('ref_name', ref_name);
    $('#id-modal-body-remove-sample').attr('tr_to_remove', tr_to_remove);
    $('#id-label-remove').val('');
  } else {
    if (attr === 'id_swap_technology_modal') {
      $('#id-label-swap').text('Do you want to swap technology for \'' + ref_name + '\'?');
      $('#id-modal-body-swap-technology').attr('pk', $(this).attr('pk'));
      $('#id-modal-body-swap-technology').attr('ref_name', ref_name);
      $('#id-modal-body-swap-technology').attr('tr_to_swap', tr_to_remove);
      $('#id-label-swap').val('');
    } else {
      if (attr === 'id_remove_unattached_modal') {
        $('#id-label-remove-unattached').text('Do you want to remove unattached samples?');
        $('#id-label-remove-unattached').val('');
      }
    }
  }
});


$('#id-remove-button').on('click', function () {

  $.ajax({
    url: $('#id-modal-body-remove-sample').attr("remove-single-value-url"),
    data: {
      sample_id: $('#id-modal-body-remove-sample').attr('pk'),
      csrfmiddlewaretoken: '{{ csrf_token }}'
    }, // data sent with the post request

    success: function (data) {

      if (data['is_ok']) {

        /// remove line
        document.getElementById($('#id-modal-body-remove-sample').attr('tr_to_remove')).remove();

        /// add message with informaton
        $('#id_messages_remove').append('<div class="alert alert-dismissible alert-success">' +
          'The sample \'' + $('#id-modal-body-remove-sample').attr('ref_name') + '\' was successfully removed.' +
          '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
          '</div>');

        /// remove one for id id-total-list
        var total_samples = $('#id-total-list').text().split(" ");
        if (total_samples.length == 3) {
          $('#id-total-list').text('Total samples: ' + (parseInt(total_samples[2]) - 1))
        }
      }
      else {
        /// special message in case present_in_televir_project == True
        if (data['present_in_televir_project']) {
          $('#id_messages_remove').append('<div class="alert alert-dismissible alert-warning">' +
            'The sample \'' + $('#id-modal-body-remove-sample').attr('ref_name') + '\' was not removed,' +
            ' because it is present in a Televir project.' +
            '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
            '</div>');
        } else {
          /// add message with informaton
          $('#id_messages_remove').append('<div class="alert alert-dismissible alert-warning">' +
            'The sample \'' + $('#id-modal-body-remove-sample').attr('ref_name') + '\' was not removed.' +
            '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
            '</div>');

        }
      }
    },

    // handle a non-successful response
    error: function (xhr, errmsg, err) {
      alert(errmsg);
      console.log(xhr.status + ": " + xhr.responseText); // provide a bit more info about the error to the console
    }
  });
});


$('#id-remove-unattached-button').on('click', function () {

  $.ajax({
    url: $('#id-modal-body-remove-unattached').attr("remove-single-value-url"),
    data: {
      csrfmiddlewaretoken: '{{ csrf_token }}'
    }, // data sent with the post request

    success: function (data) {

      if (data['is_ok']) {

        /// add message with informaton
        $('#id_messages_remove').append('<div class="alert alert-dismissible alert-success">' +
          data["message_number_samples_removed"] + ' Refresh the window to see the changes.' +
          '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
          '</div>');

      }
      else {
        /// special message in case present_in_televir_project == True
        /// add message with informaton
        $('#id_messages_remove').append('<div class="alert alert-dismissible alert-warning">' +
          'Unused samples were not properly removed.' +
          '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
          '</div>');

      }
    },

    // handle a non-successful response
    error: function (xhr, errmsg, err) {
      alert(errmsg);
      console.log(xhr.status + ": " + xhr.responseText); // provide a bit more info about the error to the console
    }
  });
});


$('#id-swap-button').on('click', function () {

  $.ajax({
    url: $('#id-modal-body-swap-technology').attr("swap-single-value-url"),
    data: {
      sample_id: $('#id-modal-body-swap-technology').attr('pk'),
      csrfmiddlewaretoken: '{{ csrf_token }}'
    }, // data sent with the post request

    success: function (data) {
      if (data['is_ok']) {

        /// update line info
        document.getElementById($('#id-modal-body-swap-technology').attr('tr_to_swap')).remove();
        ///document.getElementById($('#id-modal-body-swap-technology').attr('tr_to_swap'))[9].innerHTML = "Processing";

        /// add message with informaton
        $('#id_messages_remove').append('<div class="alert alert-dismissible alert-success">' +
          'The sample \'' + $('#id-modal-body-swap-technology').attr('ref_name') + '\' was successfully swapped.' +
          data['message'] +
          '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
          '</div>');

      }
      else {
        /// special message in case present_in_televir_project == True
        if (data['present_in_televir_project']) {
          $('#id_messages_remove').append('<div class="alert alert-dismissible alert-warning">' +
            'Tecnology for the sample \'' + $('#id-modal-body-swap-technology').attr('ref_name') + '\' was not swapped,' +
            ' because it is present in a Televir project.' +
            '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
            '</div>');
        } else {
          /// add message with informaton
          $('#id_messages_remove').append('<div class="alert alert-dismissible alert-warning">' +
            'Tecnology for the sample \'' + $('#id-modal-body-swap-technology').attr('ref_name') + '\' cannot be swapped.' +
            data['message'] +
            '<button type="button" class="close" data-dismiss="alert" aria-hidden="true">&times;</button>' +
            '</div>');

        }
      }
    },

    // handle a non-successful response
    error: function (xhr, errmsg, err) {
      alert(errmsg);
      console.log(xhr.status + ": " + xhr.responseText); // provide a bit more info about the error to the console
    }
  });
});