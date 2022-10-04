document.getElementById("sample_report").addEventListener("click", function (e) {
    if (e.target.tagName === "A" && e.target.id === "plotShow") {
        e.preventDefault();
        var row = e.target.parentNode.parentNode;
        while ((row = nextTr(row)) && !/\bparent\b/.test(row.className))
            toggle_it(row);
    }
});

function nextTr(row) {
    while ((row = row.nextSibling) && row.nodeType != 1);
    return row;
}

function toggle_it(item) {
    if (/\bopen\b/.test(item.className))
        item.className = item.className.replace(/\bopen\b/, " ");
    else
        item.className += " open";
}