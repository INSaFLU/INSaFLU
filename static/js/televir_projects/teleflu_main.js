
function showParametersHover(element) {
    var parameters = element.getAttribute('data-parameters');
    var parametersContainer = element.querySelector('.parameters-container-hover');
    parametersContainer.innerHTML = parameters;

    // Hide all other .parameters-container-hover elements
    var parentDiv = element.parentNode;
    var allParametersContainers = parentDiv.querySelectorAll('.parameters-container-hover');

    // Hide all .parameters-container-hover elements within the same parent div
    // Toggle visibility
    if (parametersContainer.style.display === 'none' || parametersContainer.style.display === '') {
        allParametersContainers.forEach(function(container) {
            container.style.display = 'none';
        });
    
        parametersContainer.style.display = 'block';
    } else {
        allParametersContainers.forEach(function(container) {
            container.style.display = 'none';
        });    
        parametersContainer.style.display = 'none';
    }


}

function showParameters(element) {
    var parameters = element.getAttribute('data-parameters');
    var parametersContainer = element.nextElementSibling;

    // If the parameters are already displayed, hide them
    if (parametersContainer.style.opacity == 1) {
        parametersContainer.style.opacity = 0;
        parametersContainer.style.maxWidth = '0';
        parametersContainer.classList.remove('visible');
    }
    // Otherwise, display the parameters
    else {
        parametersContainer.innerHTML = parameters;
        parametersContainer.style.opacity = 1;
        parametersContainer.style.maxWidth = '1000px'; // Adjust this value as needed
        parametersContainer.classList.add('visible');
    }
}


var load_teleflu_workflows = function () {
    var csrf = $('#open-modal-button').attr('csrf'); 
    var project_id = $('#open-modal-button').attr('project-id');
    $.ajax({
        url: $('#open-modal-button').attr('url-workflows'),
        type: 'GET',
        data : {
            project_id: project_id,
            csrfmiddlewaretoken: csrf,
        },
        success: function (data) {
            data.mapping_workflows.forEach(function(workflow) {
                var workflowContainer = $('<div>').addClass('workflow-container workflow-main').attr('project_id', data.teleflu_project_pk).attr('csrf', csrf);

                var workflowTitle = $('<div>').addClass('workflow-title');
                var workflowName = $('<span>').addClass('workflow-name').text('Workflow ' + workflow.node);
                workflowTitle.append(workflowName);
    
                var stepContainer = $('<div>').addClass('step-container');
                workflow.modules.forEach(function(step) {
                    var stepCircleDisplay = $('<div>').addClass('step-circle-display ' + step.available).attr('title', step.module).attr('data-parameters', step.parameters).text(step.short_name);
                    var parametersContainerHover = $('<div>').addClass('parameters-container-hover');
                    stepCircleDisplay.append(parametersContainerHover);
                    // Add onclick event
                    stepCircleDisplay.click(function() {
                        showParametersHover(this);
                    });

                    stepContainer.append(stepCircleDisplay);
                });
    
                var workflowInfo = $('<div>').addClass('workflow-info clearfix');
                var workflowMapped = $('<span>').addClass('workflow-mapped').text(workflow.samples_mapped + ' / ' + data.project_nsamples);
                var mapSamplesButton = $('<button>').attr('workflow', workflow.node).attr('workflow-id', workflow.pk).attr('type', 'button').addClass('mapSamplesButton btn btn-primary').attr('data-toggle', 'modal').attr('data-target', '#id_map_workflow_modal').text('Map Samples');
                workflowInfo.append(workflowMapped, mapSamplesButton);
    
                var mappingIgv = $('<div>').addClass('mapping-igv');

                if (workflow.samples_to_stack) {
                    var stackSamplesButton = $('<i>').attr('workflow', workflow.node).attr('workflow-id', workflow.pk).addClass('stackSamplesButton stack-deploy fa fa-flask').attr('data-toggle', 'modal').attr('data-target', '#id_workflow_analysis').attr('url', "{% url 'stack_igv_teleflu_workflow' %}");
                    mappingIgv.append(stackSamplesButton);
                } else {
                    var stackSamplesButton = $('<i>').addClass('stackSamplesButton fa fa-flask').css('color', '#b3b3b3');
                    mappingIgv.append(stackSamplesButton);
                }
    
                if (workflow.stacked_html_exists) {
                    var stackSamplesIgv = $('<a>').attr('href', workflow.stacked_html).addClass('stackSamplesIGV fa fa-eye').attr('target', '_blank').attr('title', 'Stacked IGV');
                    var stackSamplesVCF = $('<a>').attr('href', workflow.stacked_vcf).attr('rel', "nofolllow").attr('href', workflow.stacked_vcf).attr('download', 'stacked.vcf').addClass('stackSamplesVCF fa fa-download').attr('title',  "Download Stacked VCF"); 
                    mappingIgv.append(stackSamplesIgv);
                    mappingIgv.append(stackSamplesVCF);
                }
    
    
                workflowContainer.append(workflowTitle, stepContainer, workflowInfo, mappingIgv);
    
                $('#workflow-list').append(workflowContainer); // Append the new div to the body
            });

            $(".mapSamplesButton").click(function () {
                var workflow= $(this).attr('workflow');
                var workflow_id= $(this).attr('workflow-id');
                $('#id-label-map-workflow').text('Map Samples to Workflow ' + workflow + '?');
                $('#id-map-button').attr('workflow', workflow_id);            
            });
    
            $('.stack-deploy').click(function() {
                var project_id = $('.workflow-main').attr('project_id');
                var workflow_id = $(this).attr('workflow-id');
                var url = $(this).attr('url');
    
    
                $.ajax({
                    url: url,
                    type: 'POST',
                    data: {
                        project_id: project_id,
                        workflow_id: workflow_id,
                        csrfmiddlewaretoken: csrf
                    },
                    success: function(data) {
                        
                        if (data['is_ok'] === true) {
                            location.reload();
                        } else if (data['running'] === true) {
                            alert('Stacked Samples creation is running');
                        } else if (data['exists'] === true) {
                            alert('Stacked Samples already exists');
                        } else {
                            alert('Error creating stacked samples');
                        }
                    }
                });
            });

            
    


        }
    });
};



var addToProject = function (workflow, project_id) {
    var url = $('.add-to-project-button').attr('url');
    var csrf= $('.add-to-project-button').attr('csrf');

    var data = {
        'leaf_id': workflow,
        'project_id': project_id,
        'csrfmiddlewaretoken': csrf,
    };
    $.ajax({
        type: "POST",
        url: url,
        data: data,
        success: function (data) {
            if (data["is_ok"] === true) {
                alert('Workflow added to project');
                location.reload();
            } else if (data["exists"] === true) {
                alert('Workflow already exists in project');
            } else {
                alert('Error adding workflow to project');
            }
        }
    });
}


var buttons_background = function () {

    $("#id-insaflu-button").click(function () {
        var url = $(this).attr('utl');
        var project_id = $('.workflow-main').attr('project_id');
        var csrf = $('.workflow-main').attr('csrf');
        var data = {
            'project_id': project_id,
            'csrfmiddlewaretoken': csrf,
        };
        $.ajax({
            type: "POST",
            url: url,
            data: data,
            success: function (data) {
                if (data["is_ok"] === true) {
                    alert('INSaFLU Project linked');
                    location.reload();
                } else if (data["exists"] === true) {
                    alert('INSaFLU Project already exists');
                } else {
                    alert('Error linking INSaFLU Project');
                }
            }
        });
    });

        
    $(".mapSamplesButton").click(function () {
        var workflow= $(this).attr('workflow');
        var workflow_id= $(this).attr('workflow-id');
        $('#id-label-map-workflow').text('Map Samples to Workflow ' + workflow + '?');
        $('#id-map-button').attr('workflow', workflow_id);            
    });

    

    $("#id-map-button").click(function () {

        var workflow= $(this).attr('workflow');
        var url = $(this).attr('utl');
        var project_id = $('.workflow-main').attr('project_id');
        var csrf = $('.workflow-main').attr('csrf');
        var data = {
            'workflow_id': workflow,
            'project_id': project_id,
            'csrfmiddlewaretoken': csrf,
        };
        $.ajax({
            type: "POST",
            url: url,
            data: data,
            success: function (data) {
                if (data["is_ok"] === true) {
                    alert('Workflow mapped to project');
                    location.reload();
                } else if (data["is_empty"] === true) {
                    alert('No samples to map.');
                } else {
                    alert('Error mapping workflow to project');
                }
            }
        });
    });
}