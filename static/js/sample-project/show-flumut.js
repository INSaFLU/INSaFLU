$('#collapseFlumut').on('shown.bs.collapse', function () {
    show_flumut();
});

// Remove the table from the screen when collapsed
$('#collapseFlumut').on('hidden.bs.collapse', function () {
    $('#showflumut').empty();
});

// Function to show the flumut table
function show_flumut() {
    $.ajax({
        beforeSend: function() {
            $('#showflumut').empty();
            $('#loader_showflumut').show();
        },
        complete: function(){
            // $('#loader_showflumut').hide();
        },
        data: { 
            project_id: $('#showflumut').attr("project_id"),
            csrfmiddlewaretoken: '{{ csrf_token }}'
        }, // data sent with the get request
        url: $('#showflumut').attr("show-flumut-url"),
        success: function (data) {
            if (data['is_ok']) {
                // Create table
                d3.text(data['url_path_flumut'], function(datasetText) {
                    var parsedCSV = d3.tsv.parseRows(datasetText);

                    var showflumut = document.getElementById('showflumut');
                    
                    var content = '<table id="table_with_flumut_id" style="width: 100%;">\n<thead>\n';
                    var count = 0;
                    var number_cells = 0;
                    var number_cells_processed = 0;
                    parsedCSV.forEach(function(row) {
                        content += "<tr>";
                        number_cells_processed = 0;
                        row.forEach(function(cell) {
                            content += "<td>" + cell + "</td>" ;
                            number_cells_processed += 1;
                            if (count === 0){
                                number_cells += 1;
                            }
                        });
                        // Sometimes the lines don't have all fields
                        if (count != 0 && number_cells_processed != number_cells){
                            for (i = 0; i < (number_cells - number_cells_processed); i++) {
                                content += "<td></td>";
                            }
                        }
                        content += "</tr>";
                        
                        if (count === 0){
                            content += "</thead>\n<tbody>";
                        }
                        count += 1;
                    });
                    content += '</tbody></table>';
                    $('#showflumut').append(content);
                    
                    var tf = new TableFilter('table_with_flumut_id', {
                        base_path: data['static_table_filter'] + '/',
                        auto_filter: {
                            delay: 500 // milliseconds
                        },
                        alternate_rows: true,
                        rows_counter: true,
                        btn_reset: true,
                        msg_filter: 'Filtering...',
                        rows_counter: {
                            text: 'Variants: '
                        },
                        grid_layout: {
                            width: '100%'
                        },
                        col_types: [
                            'string', // Sample
                            'string', // Marker
                            'string', // Mutations
                            'string', // Effect
                            'string', // Subtype
                            'string', // Litterature

                        ],
                        col_widths: [
                            '15%',
                            '15%',
                            '15%',
                            '20%',
                            '15%',
                            '20%',

                        ],
                        col_1: 'select',
                        extensions: [{ name: 'sort' }],
                    });
                    tf.init();

                    // Add download possibility
                    content = '<a href="#" class="helpBtn" onclick="download_flumut_as_csv(\'table_with_flumut_id\');">Download as CSV</a>';
                    var element_rdiv = document.querySelector(".rdiv");
                    if (typeof(element_rdiv) != "undefined"){
                        let span = document.createElement('span');
                        span.innerHTML = content;
                        element_rdiv.appendChild(span);
                    }

                    // Hide loading icon
                    $('#loader_showflumut').hide();
                });
            } else {
                $('#loader_showflumut').hide();
                $('#showflumut').append('<div class="alert alert-warning alert-dismissable text-center"><strong>Fail</strong> to load the variants file.</div>');
            }
        },
        // Handle a non-successful response
        error: function(xhr, errmsg, err) {
            $('#loader_showflumut').hide();
            $('#showflumut').append('<div class="alert alert-warning alert-dismissable text-center"><strong>Fail</strong> to load the variants file.</div>');
            console.log(xhr.status + ": " + xhr.responseText); // Provide a bit more info about the error to the console
        }
    });
}

// Quick and simple export target #table_id into a csv
function download_flumut_as_csv(table_id, separator = ',') {
    // Do the Head, Select rows from table_id
    var rows = document.querySelectorAll('div.grd_headTblCont tr');
    // Construct csv
    var csv = [];
    for (var i = 0; i < rows.length; i++) {
        var row = [], cols = rows[i].querySelectorAll('td, th');
        for (var j = 0; j < cols.length; j++) {
            var data = cols[j].innerText.replace(/(\r\n|\n|\r)/gm, '').replace(/(\s\s)/gm, ' ');
            data = data.replace(/"/g, '""');
            row.push('"' + data + '"');
            break;
        }
        csv.push(row.join(separator));
    }
    
    // Do the body, Select rows from table_id
    rows = document.querySelectorAll('table#' + table_id + ' tr.odd, table#' + table_id + ' tr.even');
    for (var i = 0; i < rows.length; i++) {
        var row = [], cols = rows[i].querySelectorAll('td, td');
        for (var j = 0; j < cols.length; j++) {
            // Clean innertext to remove multiple spaces and jumpline (break csv)
            var data = cols[j].innerText.replace(/(\r\n|\n|\r)/gm, '').replace(/(\s\s)/gm, ' ');
            // Escape double-quote with double-double-quote
            data = data.replace(/"/g, '""');
            // Push escaped string
            row.push('"' + data + '"');
        }
        csv.push(row.join(separator));
    }
    var csv_string = csv.join('\n');
    // Download it
    var filename = 'export_' + table_id + '_' + new Date().toLocaleDateString() + '.csv';
    var link = document.createElement('a');
    link.style.display = 'none';
    link.setAttribute('target', '_blank');
    link.setAttribute('href', 'data:text/csv;charset=utf-8,' + encodeURIComponent(csv_string));
    link.setAttribute('download', filename);
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
}