

/* televir_projects/sample_detail.html
*/

/// Display inner table if it is not active, otherwise hide.

$("tr.parent").find("A#plot_show").click(function(e) {

    if (e.target.tagName === "A" && e.target.id === "plot_show") {
        e.preventDefault();
        var parent_tr = e.target.parentNode.parentNode;
        console.log(parent_tr);
        var child_tr = nextTr(parent_tr);
        console.log(child_tr);
        $(child_tr).toggleClass("active");
        
    }   

});

/// deploy remapping for single reference

$(document).on("click", "a", function (e) {
    var id= $(this).attr("id");
    var ref_id= $(this).attr("ref_id");
    if (id === "remap_reference") {
        console.log(ref_id);
        $.ajax({
            url: "{% url 'deploy_televir_map' %}",
            type: "POST",
            data: {
                'csrfmiddlewaretoken': '{{ csrf_token }}',
                'reference_id': ref_id
            },
            success: function (data) {
                console.log(data);
                if (data["is_ok"] === true) {
                    alert("Reference remap deployed");
                } else {
                    alert("Reference remapping failed");
                }
            },
            error: function (data) {
                console.log(data);
                alert("Reference remapping failed");
            }
        });
    }
});


function nextTr(row) {
while ((row = row.nextSibling) && row.nodeType != 1);
return row;
}

function toggle_it(item) {
    if (item.className != "igv_tab")
        if (/\bopen\b/.test(item.className))
            item.className = item.className.replace(/\bopen\b/, " ");
        else
            item.className += " open";
}



/// Download result table as csv (summary statistics only).

function download_table_as_csv(table_id, separator = ',') {
    // Select rows from table_id
    var rows = document.querySelectorAll('table#' + table_id + ' tr');
    // Construct csv
    var csv = [];
    for (var i = 0; i < rows.length; i++) {
        var row = [],
            cols = rows[i].querySelectorAll('td, th');
                    
        if ( rows[i].className == "header" ) {
        
            for (var j = 0; j < cols.length; j++) {
                // Clean innertext to remove multiple spaces and jumpline (break csv)
                var data = cols[j].innerText.replace(/(\r\n|\n|\r)/gm, '').replace(/(\s\s)/gm, ' ')
                // Escape double-quote with double-double-quote (see https://stackoverflow.com/questions/17808511/properly-escape-a-double-quote-in-csv)
                data = data.replace(/"/g, '""');
                // Push escaped string
                row.push(data);
            }
            csv.push(row.join(separator));
        }

        if ( rows[i].className == "parent" ) {
        
            for (var j = 0; j < cols.length; j++) {
                // Clean innertext to remove multiple spaces and jumpline (break csv)
                var data = cols[j].innerText.replace(/(\r\n|\n|\r)/gm, '').replace(/(\s\s)/gm, ' ')
                // Escape double-quote with double-double-quote (see https://stackoverflow.com/questions/17808511/properly-escape-a-double-quote-in-csv)
                data = data.replace(/"/g, '""');
                // Push escaped string
                row.push(data);
            }
            csv.push(row.join(separator));
        }
    }
    var csv_string = csv.join('\n');
    // Download it
    var filename = 'export_' + table_id + '_' + new Date().toLocaleDateString() + '.csv';
    var link = document.createElement('a');
    //link.style.display = 'none';
    link.setAttribute('target', '_blank');
    link.setAttribute('href', 'data:text/csv;charset=utf-8,' + encodeURIComponent(csv_string));
    link.setAttribute('download', filename);
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
}