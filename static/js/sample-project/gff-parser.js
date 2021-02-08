


var gff = function() {};

/// parseLines from text
gff.parseLines = function(file, last_name_sequence) {
	var lines = file.split("\n");
	var config = {};
	var arr = [];
	var line_number = 0
	var colors = ["Chocolate", "CornflowerBlue", "DarkKhaki", "DarkSeaGreen"];
	config.type = "gff3"

	for (var i = 0; i < lines.length; i++) {
		// ignore comments for now
		var line = lines[i];
		if (line.length === 0 || line[0] === "#")
			continue;

		line = gff.parseLine(line, last_name_sequence, colors, line_number);
		if (line !== undefined)
			line_number += 1;
			arr.push(line);
	}
	
	return {
		features: arr,
		config: config
	};
};


/*
 * parses one GFF line and returns it
 */
gff.extractKeys = function extractKeys(attr) {
	// extract key-value definitions
	var attributes = {};
	var attrArr = attr.split(";");
	attrArr.forEach(function(el) {
		var keyArr, key, val;
		if (el.indexOf("=") > 0) {
			keyArr = el.split("=");
			key = keyArr[0];
			val = keyArr[1];
			attributes[key] = val;
		} else if (el.indexOf(" ") > 0) {
			keyArr = el.split(" ");
			key = keyArr[0];
			val = keyArr[1].replace(/"/g, '');
			attributes[key] = val;
		}
	});
	return attributes;
};

/**
 * parses GFF and returns a dictionary of all seqs with their features
 * @method parseSeqs
 * @param {String} file GFF file
 * @return {String} Returns dictionary of sequences with an array of their features
 */
gff.parseSeqs = gff.parse = function(file, last_name_sequence) {
	var obj = gff.parseLines(file, last_name_sequence);
	var seqs = {};
	obj.features.forEach(function(entry) {
		var key = entry.seqname;
		if (seqs[key] === undefined) seqs[key] = [];
		delete entry.seqname;
		seqs[key].push(entry);
	});
	delete obj.features;
	obj.seqs = seqs;
	return obj;
};

/*
 * parses one GFF line and returns it
 */
gff.parseLine = function(line, last_name_sequence, colors, number_line) {
	var tLine = {};
	var columns = line.split(/\s+/);
	// ignore empty lines
	if (columns.length === 1)
		return;

	tLine.seqname = last_name_sequence; //columns[0];
	tLine.source = columns[1];
	tLine.feature = columns[2];
	tLine.start = parseInt(columns[3]);
	tLine.end = parseInt(columns[4]);
	tLine.score = columns[5]; 	// only DNA, RNA
	tLine.strand = columns[6]; 	// only DNA, RNA
	tLine.frame = columns[7]; 	// only DNA, RNA
	tLine.fillColor = colors[number_line % colors.length];
	tLine.text = columns[2];
	var attr = columns.slice(8).join(" "); // plain text comments

	// remove undefined (dot)
	Object.keys(tLine).forEach(function(key) {
		if (typeof(tLine[key]) === "string") {
			tLine[key] = tLine[key].trim(); // triming is important
		}
		if (tLine[key] === ".") {
			tLine[key] = undefined;
		}
	});

	// parse optional parameters
	if (tLine.score) {
		tLine.score = parseFloat(tLine.score);
	}
	if (tLine.frame) {
		tLine.frame = parseInt(tLine.frame);
	}
	tLine.attributes = gff.extractKeys(attr);
	
	if (tLine.attributes.hasOwnProperty("gene")){
		var str_gene_name = '   --   ' + tLine.attributes['gene'];
		tLine.text = str_gene_name.repeat((tLine.end - tLine.start) / str_gene_name.length )
	}
	return tLine;
};

