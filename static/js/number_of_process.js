
/// 
function get_process(){
	var message = "Jobs - To Process (0)   Running (0)   Global (0)";
	$('#info_pocess_id').text(message)
	// get if there's any checked in the server
	$.ajax({
		url: '/managing_files/ajax/get_process_running',
		data : { 
			csrfmiddlewaretoken: '{{ csrf_token }}' }, // data sent with the post request
		success: function (data) {
			if (data['is_ok']) {
				$('#info_pocess_id').text("Jobs - To Process (" + data['process_to_run'] + ")   Running (" + data['process_running'] +
						")   Global (" + data['process_to_run_total'] + ")")
			}
			else{	// clean information
				$('#info_pocess_id').text(message)
			}
		},
		// better for a slow networks
		complete:function(data){
			setTimeout(get_process, 10000);	// 10 seconds
		}
	});
}

$().ready(function(){
	get_process();
});

