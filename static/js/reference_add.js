
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

/// only for test sge in my account
$('#id-submit-sge').on('click', function(){

	$.ajax({
        url: $('#id-submit-sge').attr("remove-single-value-url"),
        data : { 
    		csrfmiddlewaretoken: '{{ csrf_token }}'
        }, // data sent with the post request
        		
        success: function (data) {
          
        },
        
        // handle a non-successful response
        error : function(xhr,errmsg,err) {
            alert(errmsg);
            console.log(xhr.status + ": " + xhr.responseText); // provide a bit more info about the error to the console
        }
	});
});