/// JS Heatmap function

function createHeatmap(data, divId) {
    data = JSON.parse(data);
    var xValues = data.map(row => row["x"]);
    var yValues = data.map(row => row["y"]);
    var zValues = data.map(row => row["value"]);

    var trace = {
        x: xValues,
        y: yValues,
        z: zValues,
        type: 'heatmap',
        hoverongaps: false,
        colorscale: 'YlOrRd',
        zmin: 0,  // Minimum color value
        zmax: 1,  // Maximum color value
    };

    var layout = {
        annotations: [],
        width: 1000,
        margin: {
            l: 100,
            r: 100,
            b: 40,
            t: 30,
            pad: 4
        },
    };

    for ( var j = 0; j < xValues.length; j++ ) {
    var currentValue = zValues[j];
    if (currentValue < 0.6) {
        var textColor = 'white';
    }else{
        var textColor = 'black';
    }
    var result = {
        xref: 'x1',
        yref: 'y1',
        x: xValues[j],
        y: yValues[j],
        text: zValues[j].toFixed(3),
        font: {
        family: 'Arial',
        size: 13,
        },
        showarrow: false,
        font: {
        color: textColor
        }
    };
    layout.annotations.push(result);
    }

    var data = [trace];
    
    Plotly.newPlot(divId, data, layout);

    var myPlot = document.getElementById(divId);
    console.log(myPlot);
    var hoverInfo = document.getElementById('hoverinfo');
    myPlot.on('plotly_hover', function (data) {
        var infotext = data.points.map(function (d) {
            return ('x=' + d.x + ', y=' + d.y + ', value=' + d.z);
        });
        
        var x_label = data.points[0].x;
        var y_label = data.points[0].y;
        var rowIdX = "row_" + x_label;

        var rowElement = document.getElementById(rowIdX);
        if (rowElement) {
            rowElement.style.backgroundColor = "lightgray";  // Change to your preferred highlight color
        }

        var rowIdY = "row_" + y_label;
        var rowElement = document.getElementById(rowIdY);
        if (rowElement) {
            rowElement.style.backgroundColor = "lightgray";  // Change to your preferred highlight color
        }

    })
        .on('plotly_unhover', function (data) {
        
            var x_label = data.points[0].x;
            var y_label = data.points[0].y;
            var rowIdX = "row_" + x_label;
            var rowElement = document.getElementById(rowIdX);
            if (rowElement) {
                rowElement.style.backgroundColor = "";  // Change to your preferred highlight color
            }
            
            var rowIdY = "row_" + y_label;
            var rowElement = document.getElementById(rowIdY);
            if (rowElement) {
                rowElement.style.backgroundColor = "";  // Change to your preferred highlight color
            }
        }
    );

};
