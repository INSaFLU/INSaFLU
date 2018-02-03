
/// 
function get_process(){
	var message = "Jobs - To Process (0)   Running (0)   Global (0)";
	
	var name = $("#id_user_name").text().replace(/^\s+|\s+$/g, '');
	var lst_split = name.split(" ");
	if ((lst_split.length == 3 & lst_split[0] == 'demo' & lst_split[1] == '-' & lst_split[2] == 'Logout') |
		(lst_split.length == 1 & lst_split[0] == '')){
		$('#info_pocess_id').text(message)
	}
	else{
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
}

$().ready(function(){
	var message = "Jobs - To Process (0)   Running (0)   Global (0)";
	$('#info_pocess_id').text(message)
	get_process();
});

