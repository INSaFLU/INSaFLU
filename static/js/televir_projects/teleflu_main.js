
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
    var url_stack = $('#open-modal-button').attr('url-stack');
    $.ajax({
        url: $('#open-modal-button').attr('url-workflows'),
        type: 'GET',
        data : {
            project_id: project_id,
            csrfmiddlewaretoken: csrf,
        },
        success: function (data) {
            data.mapping_workflows.forEach(function (workflow) {

                var workflowContainerMain = $('<div>').addClass('workflow-container-main');
                
                var workflowContainerAction = $('<div>').addClass('workflow-container workflow-main').attr('project_id', data.teleflu_project_pk).attr('csrf', csrf);
                
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
                // mapping summary: mapping fail / mapping success / total samples
                var workflowMapped = $('<span>').addClass('workflow-mapped').text(workflow.mapped_fail + ' / '  + workflow.mapped_success + ' / ' + data.project_nsamples).attr('title', 'Mapping Fail / Mapping Success / Total Samples');
                var workflowButton = $('<button>').attr('type', 'button').addClass('workflow-summary btn btn-primary').attr('node', workflow.node);
                workflowButton.append('<i class="fa fa-list-alt large-icon" aria-hidden="true"></i>'); // Add an icon to the button
                


                var mapSamplesButton = $('<button>').attr('workflow', workflow.node).attr('workflow-id', workflow.pk).attr('type', 'button').addClass('mapSamplesButton btn btn-primary').attr('data-toggle', 'modal').attr('data-target', '#id_map_workflow_modal').text('Map Samples');
                
                mapSamplesButton.prop('disabled', true);

                if (workflow.left_to_map == true) {
                    mapSamplesButton.prop('disabled', false);
                }
                
                if (workflow.running_or_queued == true) {
                    mapSamplesButton.prop('disabled', true);
                }
                
                workflowInfo.append(workflowMapped, workflowButton, mapSamplesButton);
    
                var mappingIgv = $('<div>').addClass('mapping-igv');

                if (workflow.samples_to_stack) {
                    var stackSamplesButton = $('<i>').attr('workflow', workflow.node).attr('workflow-id', workflow.pk).addClass('stackSamplesButton stack-deploy fa fa-flask').attr('data-toggle', 'modal').attr('data-target', '#id_workflow_analysis').attr('url', url_stack).attr('title', 'Stack Samples');
                    mappingIgv.append(stackSamplesButton);
                } else {
                    var stackSamplesButton = $('<i>').addClass('stackSamplesButton fa fa-flask').css('color', '#b3b3b3');
                    mappingIgv.append(stackSamplesButton);
                }
                
                if (workflow.mapped_success > 0) {
                    var workflowIgvButton = document.createElement('a');
                    workflowIgvButton.href = 'teleflu_workflow_igv/' + workflow.pk;
                    workflowIgvButton.className = 'stackSamplesIGV workflowIGVButton';
                    workflowIgvButton.target = '_blank';
                    workflowIgvButton.title = 'Workflow IGV';
                    workflowIgvButton.style.fontWeight = 'bold'; // Make the text bold
                    //workflowIgvButton.style.fontSize = 'small'; // Make the text smaller
                    workflowIgvButton.textContent = 'IGV';
                    workflowIgvButton.style.backgroundColor = 'transparent'; // Remove the background color
                    workflowIgvButton.setAttribute('aria-hidden', 'true');
                    mappingIgv.append(workflowIgvButton);

                    // button to download a zip file with all mapping files 
                    var downloadMappingButton = document.createElement('a');
                    downloadMappingButton.href = 'teleflu_mapping_files/' + workflow.pk;
                    downloadMappingButton.className = 'downloadMappingButton fa fa-download';
                    downloadMappingButton.title = 'Download Mapping Files';
                    downloadMappingButton.style.fontWeight = 'bold'; // Make the text bold
                    //downloadMappingButton.style.fontSize = 'small'; // Make the text smaller
                    // downloadMappingButton.textContent = 'Download';
                    downloadMappingButton.style.backgroundColor = 'transparent'; // Remove the background color
                    downloadMappingButton.setAttribute('aria-hidden', 'true');
                    mappingIgv.append(downloadMappingButton);
                }

                if (workflow.stacked_html_exists) {
                    var stackSamplesVariantIgv = $('<a>').attr('href', workflow.stacked_html).addClass('stackSamplesIGV fa fa-barcode').attr('target', '_blank').attr('title', 'Variants VCF IGV').attr('aria-hidden', 'true');
                    var stackSamplesVCF = $('<a>').attr('href', workflow.stacked_variants_vcf).attr('rel', "nofolllow").attr('href', workflow.stacked_variants_vcf).attr('download', 'stacked.vcf').addClass('stackSamplesVCF fa fa-download').attr('title', "Download Variants VCF"); 
                    mappingIgv.append(stackSamplesVariantIgv);
                    mappingIgv.append(stackSamplesVCF);
                }
    
                workflowContainerAction.append(workflowTitle, stepContainer, workflowInfo, mappingIgv);
                if (workflow.running_or_queued == true) {
                    var spinnerContainer = $('<div class="spinner-container"></div>');
                    // Create the spinner icon
                    var spinner = $('<i id="workflow-container-spinner" class="fa fa-spin fa-circle-o-notch"></i>');
                    // Append the spinner to the div
                    spinnerContainer.append(spinner);
                    // Append the div to the workflow container
                    workflowContainerAction.append(spinnerContainer);
                }
                
                /// Summary List
                var sampleSummary = workflow.sample_summary;
                var summaryList = document.createElement('div');
                summaryList.className = 'summary-list workfow-' + workflow.node;
                summaryList.style.display = 'none'; // Hide the list initially
                // Create table and add 'summary-table' class
                var table = document.createElement('table');
                table.classList.add('summary-table');

                // Create header row
                var headerRow = document.createElement('tr');
                var headers = ['Sample', 'Mapped', 'Success', 'Coverage', 'Windows Covered', 'Depth', 'Mapped Reads', 'Start ::', 'Mapped ::', 'Error Rate'];

                headers.forEach(function(header) {
                    var th = document.createElement('th');
                    th.innerHTML = header;
                    headerRow.appendChild(th);
                });

                table.appendChild(headerRow);

                // Create data rows
                for (var sample in sampleSummary) {
                    var row = document.createElement('tr');

                    var sampleNameCell = document.createElement('td');
                    sampleNameCell.innerHTML = sample;
                    row.appendChild(sampleNameCell);
                    
                    var mappedIndicator = sampleSummary[sample].mapped
                    //var mappedIndicator = sampleSummary[sample].mapped ? '<span style="color: green;">&#x2714;</span>' : '<span style="color: red;">&#x2718;</span>';
                    var successIndicator = sampleSummary[sample].success ? '<span style="color: green;">&#x2714;</span>' : '<span style="color: red;">&#x2718;</span>';
                    var coverageIndicator = sampleSummary[sample].coverage;
                    var windowsCoveredIndicator = sampleSummary[sample].windows_covered;
                    var depthIndicator = sampleSummary[sample].depth;
                    var mappedReadsIndicator = sampleSummary[sample].mapped_reads;
                    var start_proportion = sampleSummary[sample].start_prop;
                    var mapped_proportion = sampleSummary[sample].mapped_prop;
                    var error_rate = sampleSummary[sample].error_rate;

                    var indicatorValues = [mappedIndicator, successIndicator, coverageIndicator, windowsCoveredIndicator,
                        depthIndicator, mappedReadsIndicator, start_proportion, mapped_proportion, error_rate];

                    indicatorValues.forEach(function(value) {
                        var cell = document.createElement('td');
                        cell.innerHTML = value;
                        row.appendChild(cell);
                    });

                    table.appendChild(row);
                }

                summaryList.appendChild(table);
                
                // Create and append the download button
                var downloadBtn = document.createElement('button');
                downloadBtn.innerHTML = 'Download TSV';
                downloadBtn.className = 'download-tsv-btn';
                downloadBtn.onclick = function() {
                    var tsvData = generateTSVData(workflow.sample_summary);
                    downloadTSV('summary-table.tsv', tsvData);
                };

                // Assuming 'summaryList' is the parent element where you want to append the button
                summaryList.appendChild(downloadBtn);
                //summaryList.style.display = 'block'; // Make sure the list (and button) is visible
                
                workflowContainerMain.append(workflowContainerAction);
    
                workflowContainerMain.append(summaryList);

                $('#workflow-list').append(workflowContainerMain); // Append the new div to the body
            });


            $(".workflow-summary").click(function () {
                var node = $(this).attr('node');
                var summaryList = $('.summary-list.workfow-' + node)[0];
                // Toggle list visibility on click
                if (summaryList.style.display === 'none') {
                    summaryList.style.display = 'flex';
                } else {
                    summaryList.style.display = 'none';
                }
            });


            $(".mapSamplesButton").click(function () {
                var workflow= $(this).attr('workflow');
                var workflow_id = $(this).attr('workflow-id');
                $('#id-label-map-workflow').text('Map Samples to Workflow ' + workflow + '?');
                $('#id-map-button').attr('workflow', workflow_id);            
            });
    
            $('.stack-deploy').click(function() {
                var project_id = $('.workflow-main').attr('project_id');
                var workflow_id = $(this).attr('workflow-id');
                var url = $('#open-modal-button').attr('url-stack');
                var csrf = $('#open-modal-button').attr('csrf');

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

// Assuming this code is added at the end of your existing script

// Function to generate TSV data
function generateTSVData(sampleSummary) {
    var headers = ['Sample', 'Reference', 'Leaf', 'Mapped', 'Success', 'Coverage', 'Windows Covered', 'Depth', 'Mapped Reads', 'Start ::', 'Mapped ::', 'Error Rate'];
    var tsvContent = headers.join('\t') + '\n'; // Header row

    for (var sample in sampleSummary) {
        var rowData = [
            sample,
            sampleSummary[sample].reference,
            sampleSummary[sample].leaf,
            sampleSummary[sample].mapped ? 'Yes' : 'No',
            sampleSummary[sample].success ? 'Yes' : 'No',
            sampleSummary[sample].coverage,
            sampleSummary[sample].windows_covered,
            sampleSummary[sample].depth,
            sampleSummary[sample].mapped_reads,
            sampleSummary[sample].start_prop,
            sampleSummary[sample].mapped_prop,
            sampleSummary[sample].error_rate
        ];
        tsvContent += rowData.join('\t') + '\n'; // Add row data
    }

    return tsvContent;
}

// Function to download TSV file
function downloadTSV(filename, text) {
    var element = document.createElement('a');
    element.setAttribute('href', 'data:text/tab-separated-values;charset=utf-8,' + encodeURIComponent(text));
    element.setAttribute('download', filename);

    element.style.display = 'none';
    document.body.appendChild(element);

    element.click();

    document.body.removeChild(element);
}


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
        var project_id = $('#id-insaflu-button').attr('project-id');
        var csrf = $('#id-insaflu-button').attr('csrf');

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

    $("#id-map-button").click(function () {

        var workflow = $(this).attr('workflow');
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