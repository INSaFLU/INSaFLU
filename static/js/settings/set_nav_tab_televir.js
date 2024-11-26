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

	if (insalfu_nav_index_main_old !== index_main) {
		for (let step = 0; step < number_indexes; step++) {
			if (step == parseInt(insalfu_nav_index_main_old)) {
				/// main
				let tech_tab= $('#main_tab li:nth-child(' + (step + 1) + ') a')
				let href_main_active = tech_tab.attr('href');
				$(href_main_active).removeClass('active');
				tech_tab.removeClass('active');

				
				/// second
				let number_groups = parseInt($('#second_tab_' + step).attr('number_groups'));

				$('#second_tab_' + step).removeClass('active');

				for (let step_li = 1; step_li <= number_groups; step_li++) {
					let group_tab = $('#second_tab_' + step + ' li:nth-child(' + step_li + ') a');

					if (group_tab.hasClass('active')) {
						let index_group = group_tab.attr('index_group');
						let href_group_active = group_tab.attr('href');
						$(href_group_active).removeClass('active');
						group_tab.removeClass('active');

						let number_pipelines = parseInt($('#third_tab_' + step + '_' + index_group).attr('number_pipelines'));
						for (let step_li_third = 1; step_li_third <= number_pipelines; step_li_third++) {
							let pipeline_tab = $('#third_tab_' + step + '_' + (step_li - 1) + ' li:nth-child(' + step_li_third + ') a');

							
							if (pipeline_tab.hasClass('active')) {
								let href_active_third = pipeline_tab.attr('href');
								$(href_active_third).removeClass('active');
								pipeline_tab.removeClass('active');
								break;
							}
						}
						/// third
						break;
					}
				}
			}
		}
		insalfu_nav_index_main_old = index_main;
		index_group_old = "0";

		$('#main_tab li:nth-child(' + (parseInt(index_main) + 1) + ') a').tab('show');
		$('#second_tab_' + index_main + ' li:first-child a').tab('show');
		$('#third_tab_' + index_main + '_0 li:first-child a').tab('show');
	}
});

$('.second_tab a').on('click', function (e) {
	e.preventDefault()
	var index_tech = $(this).attr('index_tech');
	var this_id = '#second_tab_' + index_tech;
	var index_group = $(this).attr('index_group');
	var number_groups = $(this).attr('number_groups');

	if (index_group_old !== index_group) {
		for (let step = 0; step < number_groups; step++) {
			if (step == parseInt(index_group_old)) {
				/// second
				let href_group_active = $(this_id + ' li:nth-child(' + (step + 1) + ') a').attr('href');
				$(href_group_active).removeClass('active');
				$(this_id + ' li:nth-child(' + (step + 1) + ') a').removeClass('active');
				/// third
				let number_pipelines = parseInt($('#third_tab_' + index_tech + '_' + step).attr('number_pipelines'));
				for (let step_li = 1; step_li <= number_pipelines; step_li++) {
					if ($('#third_tab_' + index_tech + '_' + step + ' li:nth-child(' + step_li + ') a').hasClass('active')) {
						let href_active = $('#third_tab_' + index_tech + '_' + step + ' li:nth-child(' + step_li + ') a').attr('href');
						$(href_active).removeClass('active');
						$('#third_tab_' + index_tech + '_' + step + ' li:nth-child(' + step_li + ') a').removeClass('active');
						break;
					}
				}
				break;
			}
		}
		index_group_old = index_group;
		
		$(this_id + 'li:nth-child(' + (parseInt(index_group) + 1) + ') a').tab('show');
		$('#third_tab_' + index_tech + '_' + index_group).tab('show');
		$('#third_tab_' + index_tech + '_' + index_group + ' li:first-child a').tab('show'); // Select first tab in third
	}
});




$(document).ready(function(){
	$('#main_tab li:first-child a').tab('show') // Select first tab in main
	$('#second_tab_0 li:first-child a').tab('show') // Select first tab in second
	$('#third_tab_0_0 li:first-child a').tab('show') // Select first tab in third
	
	/// index selected in main nav
	insalfu_nav_index_main_old = "0";
	index_group_old= "0";
});	