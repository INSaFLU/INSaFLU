function hashCode(str) { // java String#hashCode
    let hash = 0;
    for (let i = 0; i < str.length; i++) {
        hash = str.charCodeAt(i) + ((hash << 5) - hash);
    }
    return hash*130;
}

function intToRGB(i){
    let c = (i & 0x00FFFFFF)
        .toString(16)
        .toUpperCase();

    return "00000".substring(0, 6 - c.length) + c;
}

function label2number(label1, label2, dateTimeline, timeCategory, metadata, id) {

    let result = {};

    for (let i = 0; i < Object.keys(metadata).length; i++) {

        if (dateTimeline === undefined || timeCategory === undefined) {

            result[metadata[i][label1] + metadata[i][label2]] = (result[metadata[i][label1] + metadata[i][label2]] + 1) || 0;

        } else {
            if (dateTimeline >= (metadata[i][timeCategory]).split('-').join('')) result[metadata[i][label1] + metadata[i][label2]] = (result[metadata[i][label1] + metadata[i][label2]] + 1) || 0;
        }
    }

    return (result[id] || 0)
}

function numberNonEmpty(label1, label2, metadata, metaaa) {

    let result = 0;

    for (let i = 0; i < Object.keys(metadata).length; i++) {

        if (metaaa === undefined) {

            if (isNaN(metadata[i][label1]) === false && metadata[i][label1] !== "" && isNaN(metadata[i][label2]) === false && metadata[i][label2] !== "") result = result + 1;

        } else {

            if (isNaN(metadata[i][label1]) === false && metadata[i][label1] !== "" && isNaN(metadata[i][label2]) === false && metadata[i][label2] !== "" && isNaN((metadata[i][metaaa]).toString().split('-').join('')) === false) result = result + 1;

        }


    }
    return result
}

function convertTime(number) {

    let renew = number.split('');

    renew.forEach(function (value, index) {

        if ((index === 3 || index === 5) && renew.length !== 4) {
            renew[index] = renew[index] + "-";
        }
    });

    return renew.join('');
}