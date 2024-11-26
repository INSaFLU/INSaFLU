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
        colorscale: [
            ['0.0', 'rgb(49,54,149)'],
            ['0.111111111111', 'rgb(69,117,180)'],
            ['0.222222222222', 'rgb(116,173,209)'],
            ['0.333333333333', 'rgb(171,217,233)'],
            ['0.444444444444', 'rgb(224,243,248)'],
            ['0.555555555556', 'rgb(254,224,144)'],
            ['0.666666666667', 'rgb(253,174,97)'],
            ['0.777777777778', 'rgb(244,109,67)'],
            ['0.888888888889', 'rgb(215,48,39)'],
            ['1.0', 'rgb(165,0,38)']
        ],
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
    if (currentValue < 0.2 || currentValue > 0.8) {
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
    myPlot.on('plotly_hover', function (data) {

        
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



function createHeatmapClades(data, divId) {
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
        colorscale: [
            ['0.0', 'rgb(49,54,149)'],
            ['0.111111111111', 'rgb(69,117,180)'],
            ['0.222222222222', 'rgb(116,173,209)'],
            ['0.333333333333', 'rgb(171,217,233)'],
            ['0.444444444444', 'rgb(224,243,248)'],
            ['0.555555555556', 'rgb(254,224,144)'],
            ['0.666666666667', 'rgb(253,174,97)'],
            ['0.777777777778', 'rgb(244,109,67)'],
            ['0.888888888889', 'rgb(215,48,39)'],
            ['1.0', 'rgb(165,0,38)']
        ],
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
    if (currentValue < 0.2 || currentValue > 0.8) {
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
    myPlot.on('plotly_hover', function (data) {

        var x_label = data.points[0].x;
        var y_label = data.points[0].y;
        var rowIdX = "info_" + x_label;

        var rowElement = document.getElementById(rowIdX);
        if (rowElement) {
            rowElement.style.backgroundColor = "lightgray";  // Change to your preferred highlight color
        }

        var rowIdY = "info_" + y_label;
        var rowElement = document.getElementById(rowIdY);
        if (rowElement) {
            rowElement.style.backgroundColor = "lightgray";  // Change to your preferred highlight color
        }

    })
        .on('plotly_unhover', function (data) {
        
            var x_label = data.points[0].x;
            var y_label = data.points[0].y;
            var rowIdX = "info_" + x_label;
            var rowElement = document.getElementById(rowIdX);
            if (rowElement) {
                rowElement.style.backgroundColor = "";  // Change to your preferred highlight color
            }
            
            var rowIdY = "info_" + y_label;
            var rowElement = document.getElementById(rowIdY);
            if (rowElement) {
                rowElement.style.backgroundColor = "";  // Change to your preferred highlight color
            }
        }
    );

};
