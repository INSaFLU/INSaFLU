/// set wait screen
$(document).on("click", "a", function(e){
	var id_ = $(this).attr('id');
	var href = $(this).attr('href');
	var onclick = $(this).attr('onclick');
	var class_ = $(this).attr('class');
	/// if undefined set empty
	if (typeof class_ == typeof undefined){
		class_ = "";
	}
	
	// test conditions
	if (href !== '#id_set_default_modal' && onclick !== 'return false;' && id_ !== 'sidenavToggler' && !class_.startsWith("nav-") &&
		href !== 'id_turn_software_on_off'){
		wait_screen();
	}
});


$('#main_tab a').on('click', function (e) {
	e.preventDefault()
	var index_main = $(this).attr('index_main');
	var number_indexes = $(this).attr('number_indexes');

	if (insalfu_nav_index_main_old !== index_main){
		for (let step = 0; step < number_indexes; step++){
			if (step == parseInt(insalfu_nav_index_main_old)){
				/// main
				let href_main_active = $('#main_tab li:nth-child(' + (step + 1) + ') a').attr('href');
				$(href_main_active).removeClass('active');
				$('#main_tab li:nth-child(' + (step + 1) + ') a').removeClass('active');
						
				/// second
				let number_li = parseInt($('#second_tab_' + step).attr('number_pipelines'));
				for (let step_li = 1; step_li <= number_li; step_li++){
					if ($('#second_tab_' + step + ' li:nth-child(' + step_li + ') a').hasClass('active')){
						let href_active = $('#second_tab_' + step + ' li:nth-child(' + step_li + ') a').attr('href');
						$(href_active).removeClass('active');
						$('#second_tab_' + step + ' li:nth-child(' + step_li + ') a').removeClass('active');
						break;
					}
				}
				break;
			}
		}
		insalfu_nav_index_main_old = index_main;
		$('#main_tab li:nth-child(' + (parseInt(index_main) + 1) + ') a').tab('show')
		$('#second_tab_' + index_main + ' li:first-child a').tab('show')
	}
})


$(document).ready(function(){
	$('#main_tab li:first-child a').tab('show') // Select first tab in main
	$('#second_tab_0 li:first-child a').tab('show') // Select first tab in second
	
	/// index selected in main nav
	insalfu_nav_index_main_old = "0";
});	