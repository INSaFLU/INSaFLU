

/// Project creation warn if project name exists

$("#project_name").keyup(function(){
    var projectname=$(this).val();
    if(projectname!=""){
      $.ajax({
        url:'{% url 'check_project_name_exist' %}',
        type:'GET',
        data:{projectname:projectname}

      })
      .done(function(response){
        console.log(response);
        if(response=="True"){
            alert("Username already taken");    }
      })
      .fail(function(){
        console.log("failed");
      })
    }
    else{
      $(".username_error").remove();
    }
});