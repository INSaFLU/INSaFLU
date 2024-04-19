$(document).ready(function() {

    $('#id_go_back_button').on('click', function () {
        wait_screen();
    });
    
    $('.collapsible-header').click(function() {
        // precent default
        event.preventDefault();
        $('.collapsible-body').slideToggle('slow');
    });       
    
    $('.summary-header').click(function() {
        // precent default
        event.preventDefault();
        $('.collapsible-summary').slideToggle('slow');
    });      

});


/// Map Reference on click

$(document).on("click", "a", function (e) {

    if (id === "remap_reference") {
        e.preventDefault();
        var id= $(this).attr("id");
        var ref_id= $(this).attr("ref_id");
        var project_id= $(this).attr("project_id");
        var csrf_token= $("#headingMapExtra").attr("csrf-token");
        var url_remap = $("#headingMapExtra").attr("url-remap");
        

        $.ajax({
            url: url_remap,
            type: "POST",
            data: {
                'csrfmiddlewaretoken': csrf_token,
                'reference_id': ref_id,
                'project_id': project_id
            },
            success: function (data) {
                if (data["is_ok"] === true) {
                    alert("Reference remap deployed");
                } else {
                    alert("Reference remapping failed");
                }
            },
            error: function (data) {
                alert("Reference remapping failed");
            }
        });
    }
});

/// Show Mapping Subsections

function nextTr(row) {
    while ((row = row.nextSibling) && row.nodeType != 1);
    return row;
}



$("tr.parent").find("A#plot_show").click(function(e) {

    if (e.target.tagName === "A" && e.target.id === "plot_show") {
        e.preventDefault();
        var parent_tr = e.target.parentNode.parentNode;
        var child_tr = nextTr(parent_tr);
        $(child_tr).toggleClass("active");
    }
  });


  /// show heatmap
  $("a").click(function(event) {

    if (event.target.id == "showButton") {
        event.preventDefault(); 
        var clickedLink = $(this);
        var rowElement = clickedLink.closest("tr");
        
        // Access the row element or perform further actions
        // Get the next row element
        var nextRowElement = rowElement.next();
        var nextRowElementClass = nextRowElement.attr('class');
        if (nextRowElementClass == "row-to-toggle") {
            // toggle slow
            nextRowElement.toggle('slow');
        }
    }

  });


/// IGV Display functions and actions

$(".igv_browse").on('click', function (e) {
    var accid = $(this).attr('accid');
    var igv_display = document.getElementById('igv_display_' + accid);
    var igv_display_className = igv_display.className;
    
    if (/\bopen\b/.test(igv_display_className)) {
        igv_display.className = igv_display.className.replace(" open", '');
        setTimeout(function () {
            var show_igv_div = document.getElementById('show_igv_' + accid);
            show_igv_div.innerHTML = "";
        }, 300);

    } else {
        igv_display.className += " open";
        show_igv($(this));
    }
});

function replace_igv_div(accid) {
    var show_igv_div = document.getElementById('show_igv_' + accid);
    
    show_igv_div.innerHTML = "";
}


function show_igv(item) {
    var accid=item.attr('accid');
    var project_pk= item.attr('project_pk');
    var sample_pk= item.attr('sample_pk');
    var run_pk = item.attr('run_pk');
    var unique_id = item.attr('reference_id');
    var url = item.attr('show-igv-url');

    $.ajax({
        /// spin 
        beforeSend: function() {
            $('#igv_display_' + accid).show();
        },
        complete: function(){
            $('#igv_display_' + accid).hide();
        },
        
        data : { 
            'project_pk': project_pk,
            'sample_pk': sample_pk,
            'run_pk': run_pk,
            'unique_id': unique_id,
            'accid': accid,
        }, // data sent with the get request
        
        url: url,

        success: function (data) {
            if (data['is_ok']) {
                
                /// set the files names
                $('#bam_file_id').empty()
                $('#bam_file_id').append(data['bam_file_id'])
                $('#bai_file_id').empty()
                $('#bai_file_id').append(data['bai_file_id'])
                $('#vcf_file_id').empty()
                $('#vcf_file_id').append(data['vcf_file_id'])
                $('#reference_id').empty()
                $('#reference_id').append(data['reference_id'])
                $('#reference_index_id').empty()
                $('#reference_index_id').append(data['reference_index_id'])
                
                /// set options 
                options = {
                        showNavigation: true,
                        showRuler: true,
                        showChromosomeWidget: true,
                        reference: {
                        id: data['reference_name'],
                        fastaURL: data['path_reference'],
                        indexURL: data['path_reference_index'],
                    },
                    trackDefaults: {
                        bam: {
                            coverageThreshold: 0.2,
                            coverageQualityWeight: true
                        }
                    },
                        tracks: [
                        {
                            type: "bam",
                            url: data['path_bam'],
                            height: 500,
                            autoHeight: false,
                            viewAsPairs: false,
                            name: data['sample_name'],
                            colorBy: "firstOfPairStrand",
                        },                        
                    ]
                }

                browser = igv.createBrowser(document.getElementById('show_igv_' + accid), options);
            }
            else{
                $('#show_igv_' + accid).append('<div class="alert alert-warning alert-dismissable"><strong>Fail</strong> to load bam file.</div>')
            }
        },
        
        // handle a non-successful response
        error : function(xhr,errmsg,err) {
            $('#show_igv_' + accid).append('<div class="alert alert-warning alert-dismissable"><strong>Fail</strong> to load bam file.</div>')
            console.log(xhr.status + ": " + xhr.responseText); // provide a bit more info about the error to the console
        }
    });
}



//  Download the table as a CSV file

document.getElementById("report_download").addEventListener("click", function (e) {
    e.preventDefault();
    download_table_as_csv("report_table")
})

function download_table_as_csv(table_id, separator = '\t') {
    // Select rows from table_id
    var rows = document.querySelectorAll('table#' + table_id + ' tr');
    var sample_id = $("#igv_browse").attr("sample_name");
    // Construct csv
    var csv = [];
    for (var i = 0; i < rows.length; i++) {
        var row = [],
            cols = rows[i].querySelectorAll('td, th');
                    
        if ( rows[i].className == "header" ) {
        
            for (var j = 0; j < cols.length; j++) {
                // Clean innertext to remove multiple spaces and jumpline (break csv)
                var data = cols[j].innerText.replace(/(\r\n|\n|\r)/gm, '').replace(/(\s\s)/gm, ' ');

                // Escape double-quote with double-double-quote (see https://stackoverflow.com/questions/17808511/properly-escape-a-double-quote-in-csv)
                data = data.replace(/"/g, '""');
                // Push escaped string
                row.push(data);
            }
            row.push("Group");
            row.push("Best_Cov");
            csv.push(row.join(separator));
        }

        if (rows[i].classList.contains('parent')) {

            var groupName = rows[i].getAttribute('group-name');
            var first = rows[i].getAttribute('first');
        
            for (var j = 0; j < cols.length; j++) {
                // Clean innertext to remove multiple spaces and jumpline (break csv)
                var data = cols[j].innerText.replace(/(\r\n|\n|\r)/gm, '').replace(/(\s\s)/gm, ' ')
                // Escape double-quote with double-double-quote (see https://stackoverflow.com/questions/17808511/properly-escape-a-double-quote-in-csv)
                data = data.replace(/"/g, '""');
                // Push escaped string
                row.push(data);
            }
            row.push(groupName);
            row.push(first);
            csv.push(row.join(separator));
        }
    }
    var csv_string = csv.join('\n');
    // Download it
    var filename =  sample_id + 'report_date_' + new Date().toLocaleDateString() + '.tsv';
    var link = document.createElement('a');
    //link.style.display = 'none';
    link.setAttribute('target', '_blank');
    link.setAttribute('href', 'data:text/csv;charset=utf-8,' + encodeURIComponent(csv_string));
    link.setAttribute('download', filename);
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
}