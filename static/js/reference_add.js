
/// 
$().ready(function(){
	var name = $("#id_user_name").text().replace(/^\s+|\s+$/g, '');
	var lst_split = name.split(" ");
	if (lst_split.length == 3 & lst_split[0] == 'demo' & lst_split[1] == '-' & lst_split[2] == 'Logout'){
		$("#id_reference_fasta").attr("disabled", "disabled");
		$("#id_reference_genbank").attr("disabled", "disabled");
	}
	else{
		$("#id_reference_fasta").removeAttr("disabled");
		$("#id_reference_genbank").removeAttr("disabled");
	}
});

