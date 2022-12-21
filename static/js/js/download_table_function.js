function exportData(html,filename) {
    /* Get the HTML data using Element by Id */
    var csv = []
    /* Declaring array variable */
    var rows = document.querySelector('table tr');

    //iterate through rows of table
    for (var i = 0, i<rows.length; i++) {
        var row = []
        var cols = rows[i].querySelectorAll('td, th');

        for(var j = 0, j<cols.length; j++) {
            row.push(cols[j].innerText);
            csv.push(row.join(","))
        }
    }

}