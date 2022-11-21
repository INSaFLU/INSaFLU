let treeType = "circular"; // Initial shape of the phylogenetic window.PhylocanvasTree.
let color = {};
let metaNumber;
let insaphylogeo_source_images = "/static/img/insaphylogeo/"

// -------------------------------------- Metadata Dependent Operations --------------------------------------

function phylTree(metaData, data_input, root_name) {
// --- TODO set root with root_name
// -------------------------------- Generating Metadata Options Coloring Nodes

	/// call the initialization
	initialization_tree_parameters();
	
    for (let j = 0; j < Object.keys(Object.values(metaData)[0]).length; j++) {

        if (j === 0) {

            let optionMetadataa = document.createElement("OPTION");

            optionMetadataa.id = "optionMetadataa";
            optionMetadataa.innerHTML= "None";

            var element = document.getElementById("selectNodeColor");
            document.getElementById("selectNodeColor").appendChild(optionMetadataa);
        }

        let optionMetadata = document.createElement("OPTION");

        optionMetadata.id = "optionMetadata"+j;
        optionMetadata.innerHTML= Object.keys(metaData[Object.keys(metaData)[0]])[j];
        optionMetadata.value = Object.keys(metaData[Object.keys(metaData)[0]])[j];
        
        document.getElementById("selectNodeColor").appendChild(optionMetadata);
    }

// -------------------------------- Listening Changes Metadata Options Coloring Nodes

    // Listening to changes in the select element of node colors which will trigger tree update.
    document.getElementById("selectNodeColor").addEventListener("change", function () {

        for (let j = 0; j < Object.keys(Object.values(metaData)[0]).length; j++) {

            if (document.getElementById("optionMetadata"+j).selected) {

                for (let i = 0; i < (window.PhylocanvasTree.leaves).length; i++) {
					if (! Object.keys(metaData).includes(window.PhylocanvasTree.leaves[i].label)){
						continue;	
					}
                	window.PhylocanvasTree.leaves[i].setDisplay({

                        leafStyle: {
                            fillStyle: color[Object.keys(metaData[window.PhylocanvasTree.leaves[i].label])[j] + 
								Object.values(metaData[window.PhylocanvasTree.leaves[i].label])[j]]
                        }
                    });
                }
            }
        }

        if (document.getElementById("optionMetadataa").selected) {

            for (let i = 0; i < (window.PhylocanvasTree.leaves).length; i++) {

                window.PhylocanvasTree.leaves[i].setDisplay({
                    //colour: 'red',
                    //shape: 'circle', // or square, triangle, star
                    //size: 3, // ratio of the base node size
                    leafStyle: {
                        //strokeStyle: '#0000ff',
                        fillStyle: "black"
                        //lineWidth: 2,
                    }
                });
            }
        }

        window.PhylocanvasTree.setNodeSize(document.getElementById("nodeSize").value);
        window.PhylocanvasTree.setTextSize(document.getElementById("textSize").value);
    });

    // -------------------------------- Generating Metadata Options Bar Tree Leafs

    for (let j = 0; j < Object.keys(Object.values(metaData)[0]).length; j++) {

        let metadataSwitchLabel = document.createElement("LABEL"); // Metadata label.
        let metadataSwitch = document.createElement("INPUT"); // Checkbox input.

        metadataSwitchLabel.id = "metadataSwitchLabel"+j;
        metadataSwitchLabel.className = "metadataSwitchLabel";
        metadataSwitchLabel.innerHTML = Object.keys(metaData[Object.keys(metaData)[0]])[j];
        metadataSwitchLabel.style.display = "block";
        metadataSwitchLabel.style.paddingLeft = "10px";
        metadataSwitchLabel.style.paddingRight = "10px";
        metadataSwitchLabel.style.cursor = "pointer";

        metadataSwitch.id = "metadataSwitch"+j;
        metadataSwitch.className = "metadataSwitch";
        metadataSwitch.type = "checkbox";
        metadataSwitch.style.cursor = "pointer";
        metadataSwitch.style.marginLeft = "10px";

        document.getElementById("treeMetadataDivInside").appendChild(metadataSwitchLabel);
        document.getElementById("metadataSwitchLabel"+j).appendChild(metadataSwitch);
    }

    // -------------------------------- Creating Meta ({metadata: values, ...})

    let meta = {};

    for (let i = 0; i < Object.keys(Object.values(metaData)[0]).length; i++) { // Picking metadata type.

        for (let j = 0; j < Object.keys(metaData).length; j++) { // Picking strain by name.

            meta[Object.keys(Object.values(metaData)[0])[i]] = meta[Object.keys(Object.values(metaData)[0])[i]] || {};
            if (Object.values(metaData[Object.keys(metaData)[j]])[i] === "" || Object.values(metaData[Object.keys(metaData)[j]])[i] === undefined || Object.values(metaData[Object.keys(metaData)[j]])[i] === " ") {
                meta[Object.keys(Object.values(metaData)[0])[i]]["No Data"] = true;
            } else {
                meta[Object.keys(Object.values(metaData)[0])[i]][Object.values(metaData[Object.keys(metaData)[j]])[i]] = true;
            }
        }
    }

    // -------------------------------- Generating Legend Dropdown

    let legendSwitchLabel = document.createElement("DIV"); // Metadata label.
    let triangle = document.createElement("IMG"); // Metadata label.
    let block = document.createElement("IMG"); // Metadata label.
    let collapse = document.createElement("IMG"); // Metadata label.

    legendSwitchLabel.id = "legendSwitchLabel" + "test";
    legendSwitchLabel.style.display = "block";
    legendSwitchLabel.style.paddingLeft = "10px";
    legendSwitchLabel.style.paddingRight = "10px";
    legendSwitchLabel.style.backgroundColor = "white";
    legendSwitchLabel.style.color = "black";
    legendSwitchLabel.style.fontSize = "14px";
    legendSwitchLabel.style.textAlign = "center";

    triangle.id = "triangle";
    triangle.className = "insaphylo_style_buttons_top_menu";
    triangle.src = insaphylogeo_source_images + "reset.png";
    collapse.id = "collapse";
    collapse.className = "insaphylo_style_buttons_top_menu";
    collapse.src = insaphylogeo_source_images + "triangle_down.png";
    block.id = "block";
    block.className = "insaphylo_style_buttons_top_menu";
    block.src = insaphylogeo_source_images + "non_block.png";

    document.getElementById("treeLegendDivInside").appendChild(legendSwitchLabel);
    document.getElementById("legendSwitchLabel" + "test").appendChild(triangle);
    document.getElementById("legendSwitchLabel" + "test").appendChild(collapse);
    document.getElementById("legendSwitchLabel" + "test").appendChild(block);

        // Metadata Categories
        for (let i = 0; i < Object.keys(Object.values(metaData)[0]).length; i++) { // Picking metadata type.

            let legendSwitchLabel = document.createElement("DIV"); // Metadata label.
            let triangle = document.createElement("IMG"); // Metadata label.

            legendSwitchLabel.id = "legendSwitchLabel" + i;
            legendSwitchLabel.className = "legendSwitchLabel";
            legendSwitchLabel.innerHTML = Object.keys(metaData[Object.keys(metaData)[0]])[i];
            legendSwitchLabel.style.display = "block";
            legendSwitchLabel.style.paddingLeft = "10px";
            legendSwitchLabel.style.paddingRight = "10px";
            legendSwitchLabel.style.opacity = "1";
            legendSwitchLabel.style.backgroundColor = "white";
            legendSwitchLabel.style.color = "black";
            legendSwitchLabel.style.fontSize = "14px";
            legendSwitchLabel.style.cursor = "pointer";

            triangle.id = "triangle" + i;
            triangle.style.display = "inline-block";
            triangle.style.paddingLeft = "10px";
            triangle.style.paddingRight = "10px";
            triangle.style.opacity = "1";
            triangle.style.cursor = "pointer";
            triangle.style.width = "10px";
            triangle.src = insaphylogeo_source_images + "triangle_down.png";

            document.getElementById("treeLegendDivInside").appendChild(legendSwitchLabel);
            document.getElementById("legendSwitchLabel" + i).appendChild(triangle);

            // Metadata Category Values
            for (let j = 0; j < Object.keys(meta[Object.keys(meta)[i]]).length; j++) {

                let legendSwitchLabelCollapsible = document.createElement("LABEL"); // Metadata label.
                let colorPicker = document.createElement("INPUT");

                legendSwitchLabelCollapsible.className = "legendSwitchLabelCollapsible"+i;
                legendSwitchLabelCollapsible.id = "legendSwitchLabelCollapsible"+i+"-"+j;
                legendSwitchLabelCollapsible.innerHTML = Object.keys(meta[Object.keys(meta)[i]])[j];
                legendSwitchLabelCollapsible.style.display = "none";
                legendSwitchLabelCollapsible.style.paddingLeft = "10px";
                legendSwitchLabelCollapsible.style.paddingRight = "10px";
                legendSwitchLabelCollapsible.style.opacity = "1";
                legendSwitchLabelCollapsible.style.cursor = "pointer";
                if (Object.keys(meta[Object.keys(meta)[i]])[j] === "No Data") {
                    legendSwitchLabelCollapsible.style.backgroundColor = "white";
                    legendSwitchLabelCollapsible.style.color = "black";
                } else {
                    legendSwitchLabelCollapsible.style.backgroundColor = "#" + intToRGB(hashCode(Object.keys(meta[Object.keys(meta)[i]])[j] + Object.keys(meta[Object.keys(meta)[i]])[j]));
                    legendSwitchLabelCollapsible.style.color = "white";
                }

                colorPicker.className = "colorPicker";
                colorPicker.id = "colorPicker"+i+"-"+j;
                colorPicker.type = "color";
                if (Object.keys(meta[Object.keys(meta)[i]])[j] === "No Data") {
                    colorPicker.value = "#FFFFFF";
                } else {
                    colorPicker.value = "#" + intToRGB(hashCode(Object.keys(meta[Object.keys(meta)[i]])[j] + Object.keys(meta[Object.keys(meta)[i]])[j]));
                }
                colorPicker.style.display = "none";
                colorPicker.style.textAlign = "center";
                colorPicker.style.width = "30px";
                colorPicker.style.height = "15px";
                colorPicker.style.color = "white";
                colorPicker.style.border = "solid";
                colorPicker.style.borderWidth = "1px";
                colorPicker.style.borderRadius = "10px";
                colorPicker.style.borderColor = "transparent";
                colorPicker.style.backgroundColor = "black";
                colorPicker.style.marginLeft = "10px";
                colorPicker.style.cursor = "pointer";

                document.getElementById("treeLegendDivInside").appendChild(legendSwitchLabelCollapsible); //legendDiv
                document.getElementById("legendSwitchLabelCollapsible"+i+"-"+j).appendChild(colorPicker);
            }
        }

    // -------------------------------- Generating Initial Color Set

    for (let i = 0; i < Object.keys(Object.values(metaData)[0]).length; i++) {

        for (let j = 0; j < Object.keys(meta[Object.keys(meta)[i]]).length; j++) {

            document.getElementById("colorPicker"+i+"-"+j).addEventListener("click", function () {

                document.getElementById("treeLegendDivInside").style.display = "block";

                document.getElementById("treeLegendDivInside").addEventListener("mouseout", function () {

                    treeLegendDivInside.style.display = "block";
                });
            })
        }
    }

    // -------------------------------- Listening Color Input

    for (let i = 0; i < Object.keys(Object.values(metaData)[0]).length; i++) {

        for (let j = 0; j < Object.keys(meta[Object.keys(meta)[i]]).length; j++) {

            document.getElementById("colorPicker"+i+"-"+j).addEventListener("change", function () {

                document.getElementById("legendSwitchLabelCollapsible"+i+"-"+j).style.backgroundColor = document.getElementById("colorPicker"+i+"-"+j).value;
                document.getElementById("legendSwitchLabelCollapsible"+i+"-"+j).style.color = "white";

                for (let ij = 0; ij < (window.PhylocanvasTree.leaves).length; ij++) { // Iterates along all the strains.
					if (! Object.keys(metaData).includes(window.PhylocanvasTree.leaves[ij].label)){
						continue;	
					}
	                if (Object.values(metaData[window.PhylocanvasTree.leaves[ij].label])[i] === document.getElementById("legendSwitchLabelCollapsible" + i + "-" + j).innerText || (Object.values(metaData[window.PhylocanvasTree.leaves[ij].label])[i] === "" && document.getElementById("legendSwitchLabelCollapsible" + i + "-" + j).innerText === "No Data")) {
	
	                    color[Object.keys(metaData[window.PhylocanvasTree.leaves[ij].label])[i] + Object.values(metaData[window.PhylocanvasTree.leaves[ij].label])[i]] = document.getElementById("colorPicker" + i + "-" + j).value;
	
	                }
	
	                if (document.getElementById("selectNodeColor").value === document.getElementById("optionMetadata" + i).innerText) {
	
	                    window.PhylocanvasTree.leaves[ij].setDisplay({
	
	                        leafStyle: {
	                            fillStyle: color[Object.keys(metaData[window.PhylocanvasTree.leaves[ij].label])[i] + Object.values(metaData[window.PhylocanvasTree.leaves[ij].label])[i]] // 2nd option input color value
	                        }
	                    });
	                }

                    metaNumber = i;

                    if (document.getElementById("metadataSwitch"+i).checked === true) {

                        (window.PhylocanvasTree.leaves[ij].data)[Object.keys(metaData[window.PhylocanvasTree.leaves[ij].label])[i]] = (window.PhylocanvasTree.leaves[ij].data)[Object.keys(metaData[window.PhylocanvasTree.leaves[ij].label])[i]] || {};

                        (window.PhylocanvasTree.leaves[ij].data)[Object.keys(metaData[window.PhylocanvasTree.leaves[ij].label])[i]]["colour"] = color[Object.keys(metaData[window.PhylocanvasTree.leaves[ij].label])[i] + Object.values(metaData[window.PhylocanvasTree.leaves[ij].label])[i]];

                        if (document.getElementById("metadataSwitchDisplay").checked === true) {

                            (window.PhylocanvasTree.leaves[ij].data)[Object.keys(metaData[window.PhylocanvasTree.leaves[ij].label])[i]]["label"] = Object.values(metaData[window.PhylocanvasTree.leaves[ij].label])[i] || "No Data";

                        } else {

                            (window.PhylocanvasTree.leaves[ij].data)[Object.keys(metaData[window.PhylocanvasTree.leaves[ij].label])[i]]["label"] = "";

                        }
                    }
                }

                    //document.getElementById("selectNodeColor").value = document.getElementById("optionMetadata" + i).innerText;

                    window.PhylocanvasTree.setNodeSize(document.getElementById("nodeSize").value);
                    window.PhylocanvasTree.setTextSize(document.getElementById("textSize").value);

                document.getElementById("treeLegendDivInside").addEventListener("mouseout", function () {

                    treeLegendDivInside.style.display = "none";
                });

                document.getElementById("triangle").style.transform = "rotate(0deg)";

            })
            }
    }

    // -------------------------------- Listening Color Reset

            document.getElementById("triangle").addEventListener("click", function () {

                    for (let j = 0; j < Object.keys(Object.values(metaData)[0]).length; j++) {

                        // Color Array Update
                        for (let i = 0; i < Object.keys(metaData).length; i++) {

                            if (Object.values(metaData[window.PhylocanvasTree.leaves[i].label])[j] === "" || Object.values(metaData[window.PhylocanvasTree.leaves[i].label])[j] === " " || Object.values(metaData[window.PhylocanvasTree.leaves[i].label])[j] === undefined) {
                                color[Object.keys(metaData[window.PhylocanvasTree.leaves[i].label])[j] + Object.values(metaData[window.PhylocanvasTree.leaves[i].label])[j]] = "white";
                            } else {
                                color[Object.keys(metaData[window.PhylocanvasTree.leaves[i].label])[j] + Object.values(metaData[window.PhylocanvasTree.leaves[i].label])[j]] = "#" + intToRGB(hashCode(Object.values(metaData[window.PhylocanvasTree.leaves[i].label])[j] + Object.values(metaData[window.PhylocanvasTree.leaves[i].label])[j]));
                            }
                        }

                        // Color Metadata Blocks Update
                        if (document.getElementById("metadataSwitch"+j).checked === true) {

                            for (let ii = 0; ii < (window.PhylocanvasTree.leaves).length; ii++) {

								if (! Object.keys(metaData).includes(window.PhylocanvasTree.leaves[ii].label)){
									continue;	
								}
                                (window.PhylocanvasTree.leaves[ii].data)[Object.keys(metaData[window.PhylocanvasTree.leaves[ii].label])[j]] = (window.PhylocanvasTree.leaves[ii].data)[Object.keys(metaData[window.PhylocanvasTree.leaves[ii].label])[j]] || {};

                                (window.PhylocanvasTree.leaves[ii].data)[Object.keys(metaData[window.PhylocanvasTree.leaves[ii].label])[j]]["colour"] = color[Object.keys(metaData[window.PhylocanvasTree.leaves[ii].label])[j] + Object.values(metaData[window.PhylocanvasTree.leaves[ii].label])[j]];

                                if (document.getElementById("metadataSwitchDisplay").checked === true) {

                                    (window.PhylocanvasTree.leaves[ii].data)[Object.keys(metaData[window.PhylocanvasTree.leaves[ii].label])[j]]["label"] = Object.values(metaData[window.PhylocanvasTree.leaves[ii].label])[j] || "No Data";

                                } else {

                                    (window.PhylocanvasTree.leaves[ii].data)[Object.keys(metaData[window.PhylocanvasTree.leaves[ii].label])[j]]["label"] = "";

                                }
                            }
                        }

                        document.getElementById("triangle").style.transform = "rotate(-30deg)";

                        // Metadata Category Values
                        for (let jj = 0; jj < Object.keys(meta[Object.keys(meta)[j]]).length; jj++) {

                            if (Object.keys(meta[Object.keys(meta)[j]])[jj] === "No Data") {
                                document.getElementById("legendSwitchLabelCollapsible"+j+"-"+jj).style.backgroundColor = "white";
                                document.getElementById("colorPicker"+j+"-"+jj).value = "#FFFFFF";
                                document.getElementById("legendSwitchLabelCollapsible"+j+"-"+jj).style.color = "black";
                            } else {
                                document.getElementById("legendSwitchLabelCollapsible"+j+"-"+jj).style.backgroundColor = "#" + intToRGB(hashCode(Object.keys(meta[Object.keys(meta)[j]])[jj] + Object.keys(meta[Object.keys(meta)[j]])[jj]));
                                document.getElementById("colorPicker"+j+"-"+jj).value = "#" + intToRGB(hashCode(Object.keys(meta[Object.keys(meta)[j]])[jj] + Object.keys(meta[Object.keys(meta)[j]])[jj]));
                            }
                        }
                    }

                if (document.getElementById("optionMetadataa").selected) {

                    for (let i = 0; i < window.PhylocanvasTree.leaves.length; i++) {

                        window.PhylocanvasTree.leaves[i].setDisplay({
                            //colour: 'red',
                            //shape: 'circle', // or square, triangle, star
                            //size: 3, // ratio of the base node size
                            leafStyle: {
                                //strokeStyle: '#0000ff',
                                fillStyle: "black"
                                //lineWidth: 2,
                            }
                        });
                    }
                } else {

                    for (let ij = 0; ij < window.PhylocanvasTree.leaves.length; ij++) { // Iterates along all the strains.

                        window.PhylocanvasTree.leaves[ij].setDisplay({

                            leafStyle: {
                                fillStyle: color[Object.keys(metaData[window.PhylocanvasTree.leaves[ij].label])[metaNumber] + Object.values(metaData[window.PhylocanvasTree.leaves[ij].label])[metaNumber]] // 2nd option input color value
                            }
                        });
                    }
                }
            });

    // -------------------------------- Listening Style Reset

    document.getElementById("resetStyleButton").addEventListener("click", function () {

        document.getElementById("nodeSize").value = "10";
        document.getElementById("textSize").value = "12";
        document.getElementById("lineWidth").value = "20";

        document.getElementById("nodeSizeLabel").innerHTML = "Node Size: " + "10" + "px";
        document.getElementById("textSizeLabel").innerHTML = "Label Size: " + "12" + "px";
        document.getElementById("lineWidthLabel").innerHTML = "Line Width: " + "1" + "px";

        window.PhylocanvasTree.lineWidth = "1";
        window.PhylocanvasTree.setNodeSize("10"); // Setting the node size accordingly to the value of the respective range element.
        window.PhylocanvasTree.setTextSize("12"); // Setting labels' font size.

        document.getElementById("resetStyleButton").style.transform = "rotate(-30deg)";
        document.getElementById("optionMetadataa").selected = true;

        for (let i = 0; i < window.PhylocanvasTree.leaves.length; i++) {

            window.PhylocanvasTree.leaves[i].setDisplay({
                //colour: 'red',
                //shape: 'circle', // or square, triangle, star
                //size: 3, // ratio of the base node size
                leafStyle: {
                    //strokeStyle: '#0000ff',
                    fillStyle: "black"
                    //lineWidth: 2,
                }
            });
        }
    });

    // -------------------------------- Listening Collapse All

    document.getElementById("collapse").addEventListener("click", function () {

        for (let i = 0; i < document.getElementsByClassName("legendSwitchLabel").length; i++) {

                for (let ii = 0; ii < document.getElementsByClassName("legendSwitchLabelCollapsible" + i).length; ii++) {

                    if (document.getElementById("collapse").src.search("img/triangle_down") !== -1) {

                        document.getElementsByClassName("legendSwitchLabelCollapsible" + i)[ii].style.display = "block";
                        document.getElementById("triangle" + i).src = insaphylogeo_source_images + "triangle_up.png";

                        if (i === document.getElementsByClassName("legendSwitchLabel").length - 1 && ii === document.getElementsByClassName("legendSwitchLabelCollapsible" + i).length - 1) {
                            document.getElementById("collapse").src = insaphylogeo_source_images + "triangle_up.png";

                        }

                    } else {

                        document.getElementsByClassName("legendSwitchLabelCollapsible" + i)[ii].style.display = "none";
                        document.getElementById("triangle" + i).src = insaphylogeo_source_images + "triangle_down.png";

                        if (i === document.getElementsByClassName("legendSwitchLabel").length - 1 && ii === document.getElementsByClassName("legendSwitchLabelCollapsible" + i).length - 1) {
                            document.getElementById("collapse").src = insaphylogeo_source_images + "triangle_down.png";
                        }
                    }
                }
        }
    });

    // -------------------------------- Listening Legend Div Block

    document.getElementById("block").addEventListener("click", function () {

        if (document.getElementById("block").src.search(insaphylogeo_source_images + "block.png") !== -1) {

            document.getElementById("block").src = insaphylogeo_source_images + "non_block.png";

            document.getElementById("treeLegendDivInside").addEventListener("mouseout", function () {

                treeLegendDivInside.style.display = "none";
            });

        } else {

            document.getElementById("block").src = insaphylogeo_source_images + "block.png";
            document.getElementById("treeLegendDivInside").addEventListener("mouseout", function () {

                treeLegendDivInside.style.display = "block";
            });

        }
    });

    // -------------------------------- Listening Style Div Block

    document.getElementById("blockStyleButton").addEventListener("click", function () {

        if (document.getElementById("blockStyleButton").src.search(insaphylogeo_source_images + "block.png") !== -1) {

            document.getElementById("blockStyleButton").src = insaphylogeo_source_images + "non_block.png";

            document.getElementById("treeStyleDivInside").addEventListener("mouseout", function () {

                treeStyleDivInside.style.display = "none";
            });

        } else {

            document.getElementById("blockStyleButton").src = insaphylogeo_source_images + "block.png";
            document.getElementById("treeStyleDivInside").addEventListener("mouseout", function () {

                treeStyleDivInside.style.display = "block";
            });

        }
    });

    // -------------------------------- Listening Clicks Expand Legend

    let aux = 0;

    for (let i = 0; i < document.getElementsByClassName("legendSwitchLabel").length; i++) {

        document.getElementById("legendSwitchLabel"+i).addEventListener("click", function () {

            for (let ii = 0; ii < document.getElementsByClassName("legendSwitchLabelCollapsible" + i).length; ii++) {

                if (document.getElementsByClassName("legendSwitchLabelCollapsible" + i)[ii].style.display === "block") {

                    document.getElementsByClassName("legendSwitchLabelCollapsible" + i)[ii].style.display = "none";
                    document.getElementById("triangle" + i).src = insaphylogeo_source_images + "triangle_down.png";
                    aux--;

                } else {

                    document.getElementsByClassName("legendSwitchLabelCollapsible" + i)[ii].style.display = "block";
                    document.getElementById("triangle" + i).src = insaphylogeo_source_images + "triangle_up.png";
                    document.getElementById("collapse").src = insaphylogeo_source_images + "triangle_up.png";
                    aux++;
                }
            }
            if(aux === 0) document.getElementById("collapse").src = insaphylogeo_source_images + "triangle_down.png";
        });
    }

    // -------------------------------- Direct Commands

    document.getElementById("phylocanvas").style.overflow = "hidden";
    document.getElementById("treeButton").src = insaphylogeo_source_images + treeType + ".png"; // Initializing the right tree button accordingly to the shape of the window.PhylocanvasTree.
    
}

// -------------------------------------- Auxiliary Function Metadata Switch --------------------------------------

function radioDetect (metaData, max) {

    // -------------------------------- Generating Initial Color Set

        for (let j = 0; j < Object.keys(Object.values(metaData)[0]).length; j++) {
        	
            for (let i = 0; i < (window.PhylocanvasTree.leaves).length; i++) {
				if (! Object.keys(metaData).includes(window.PhylocanvasTree.leaves[i].label)){
					continue;	
				}
            	if (Object.values(metaData[window.PhylocanvasTree.leaves[i].label])[j] === undefined || 
									Object.values(metaData[window.PhylocanvasTree.leaves[i].label])[j] === "") {

                    color[Object.keys(metaData[window.PhylocanvasTree.leaves[i].label])[j] + Object.values(metaData[window.PhylocanvasTree.leaves[i].label])[j]] = "white";

                } else {

                    color[Object.keys(metaData[window.PhylocanvasTree.leaves[i].label])[j] + Object.values(metaData[window.PhylocanvasTree.leaves[i].label])[j]] = "#" + intToRGB(hashCode(Object.values(metaData[window.PhylocanvasTree.leaves[i].label])[j] + Object.values(metaData[window.PhylocanvasTree.leaves[i].label])[j]));

                }
            }
        }

    // -------------------------------- Listening Changes Metadata Checkboxes

    for (let j = 0; j < max; j++) {

        document.getElementById("metadataSwitch"+j).addEventListener("change", function () {

            if (document.getElementById("metadataSwitch"+j).checked === true) {

                for (let i = 0; i < (window.PhylocanvasTree.leaves).length; i++) {
					if (! Object.keys(metaData).includes(window.PhylocanvasTree.leaves[i].label)){
						continue;	
					}
				
                    (window.PhylocanvasTree.leaves[i].data)[Object.keys(metaData[window.PhylocanvasTree.leaves[i].label])[j]] = (window.PhylocanvasTree.leaves[i].data)[Object.keys(metaData[window.PhylocanvasTree.leaves[i].label])[j]] || {};

                    (window.PhylocanvasTree.leaves[i].data)[Object.keys(metaData[window.PhylocanvasTree.leaves[i].label])[j]]["colour"] = color[Object.keys(metaData[window.PhylocanvasTree.leaves[i].label])[j] + Object.values(metaData[window.PhylocanvasTree.leaves[i].label])[j]];

                    if (document.getElementById("metadataSwitchDisplay").checked === true) {

                        (window.PhylocanvasTree.leaves[i].data)[Object.keys(metaData[window.PhylocanvasTree.leaves[i].label])[j]]["label"] = Object.values(metaData[window.PhylocanvasTree.leaves[i].label])[j] || "No Data";

                    } else {

                        (window.PhylocanvasTree.leaves[i].data)[Object.keys(metaData[window.PhylocanvasTree.leaves[i].label])[j]]["label"] = "";

                    }
                }

            } else {

                for (let i = 0; i < window.PhylocanvasTree.leaves.length; i++) {

					if (! Object.keys(metaData).includes(window.PhylocanvasTree.leaves[i].label)){
						continue;	
					}
                    delete (window.PhylocanvasTree.leaves[i].data)[Object.keys(metaData[window.PhylocanvasTree.leaves[i].label])[j]];

                }
            }

            window.PhylocanvasTree.setTreeType(treeType); // Choosing type of tree: takes radial, rectangular, circular, diagonal and hierarchy.
            window.PhylocanvasTree.setNodeSize(document.getElementById("nodeSize").value);
            window.PhylocanvasTree.setTextSize(document.getElementById("textSize").value);
        });
    }
}

function initialization_tree_parameters() {
	// -------------------------------------- Toggle Gear --------------------------------------
	
	let treeToggleGear = document.createElement("IMG");
	
	treeToggleGear.id = "toggle";
	treeToggleGear.style.position = "absolute";
	treeToggleGear.style.right = "10px";
	treeToggleGear.style.top = "10px";
	treeToggleGear.onclick = toggle;
	treeToggleGear.style.zIndex = "1";
	treeToggleGear.src = insaphylogeo_source_images + "pixelmator/gear_black.png";
	treeToggleGear.style.width = "25px";
	treeToggleGear.style.cursor = "pointer";
	
	document.getElementById("phylocanvas").appendChild(treeToggleGear);
	
	// -------------------------------------- Tree --------------------------------------
	let treeDiv = document.createElement("DIV");
	let treeButton = document.createElement("IMG");
	let treeDivInside = document.createElement("DIV");
	
	let treeType1 = document.createElement("IMG");
	let treeType2 = document.createElement("IMG");
	let treeType3 = document.createElement("IMG");
	let treeType4 = document.createElement("IMG");
	let treeType5 = document.createElement("IMG");
	
	treeDiv.id = "treeDiv";
	treeDiv.classList.add("toggle");
	treeDiv.style.position = "absolute";
	treeDiv.style.zIndex = "2000";
	treeDiv.style.right = "10px";
	treeDiv.style.bottom = "10px";
	treeDiv.style.display = "none";
	
	treeButton.id = "treeButton";
	treeButton.src = insaphylogeo_source_images + "circular.png";
	treeButton.style.width = "30px";
	treeButton.style.cursor = "pointer";
	
	treeDivInside.id = "treeDivInside";
	treeDivInside.class = "dropdown-content";
	treeDivInside.style.display = "none";
	treeDivInside.style.position = "absolute";
	treeDivInside.style.bottom = "30px";
	treeDivInside.style.right = "0px";
	treeDivInside.style.boxShadow = "0 8px 16px 0 rgba(0,0,0,0.9)";
	treeDivInside.style.backgroundColor = "black";
	treeDivInside.style.borderRadius = "10px";
	
	treeType1.id = "radial";
	treeType1.className = "insaphylo_style_buttons_left";
	treeType1.src = insaphylogeo_source_images + "radial.png";
	
	treeType2.id = "rectangular";
	treeType2.className = "insaphylo_style_buttons_left";
	treeType2.src = insaphylogeo_source_images + "rectangular.png";
	
	treeType3.id = "circular";
	treeType3.className = "insaphylo_style_buttons_left";
	treeType3.src = insaphylogeo_source_images + "circular.png";
	
	treeType4.id = "diagonal";
	treeType4.className = "insaphylo_style_buttons_left";
	treeType4.src = insaphylogeo_source_images + "diagonal.png";
	
	treeType5.id = "hierarchical";
	treeType5.className = "insaphylo_style_buttons_left";
	treeType5.src = insaphylogeo_source_images + "hierarchical.png";
	
	document.getElementById("phylocanvas").appendChild(treeDiv);
	document.getElementById("treeDiv").appendChild(treeButton);
	document.getElementById("treeDiv").appendChild(treeDivInside);
	document.getElementById("treeDivInside").appendChild(treeType1);
	document.getElementById("treeDivInside").appendChild(treeType2);
	document.getElementById("treeDivInside").appendChild(treeType3);
	document.getElementById("treeDivInside").appendChild(treeType4);
	document.getElementById("treeDivInside").appendChild(treeType5);
	
	// -------------------------------------- Style --------------------------------------
	
	const INSAPHYLOGEO_SIZE_MENU_WITH = "230px"
	const INSAPHYLOGEO_SIZE_MENU_LABEL_HEIGTH = "30px"
	
	const INSAPHYLOGEO_SIZE_FONTSIZE_MENU = "15px"
	const INSAPHYLOGEO_SIZE_FONTSIZE_ITENS_MENU = "13px"

	let treeStyleDiv = document.createElement("DIV");
	let treeStyleButton = document.createElement("BUTTON");
	let treeStyleDivInside = document.createElement("DIV");
	
	let propertiesStyleDiv = document.createElement("LABEL"); // Metadata label.
	let resetStyleButton = document.createElement("IMG"); // Metadata label.
	let blockStyleButton = document.createElement("IMG"); // Metadata label.
	
	let nodeSizeLabel = document.createElement("P");
	let nodeSize = document.createElement("INPUT");
	
	let textSizeLabel = document.createElement("P");
	let textSize = document.createElement("INPUT");
	
	let lineWidthLabel = document.createElement("P");
	let lineWidth = document.createElement("INPUT");
	
	let selectNodeColorLabel = document.createElement("P");
	let selectNodeColor = document.createElement("SELECT");
	
	treeStyleDiv.id = "treeStyleDiv";
	treeStyleDiv.classList.add("toggle");
	treeStyleDiv.style.position = "absolute";
	treeStyleDiv.style.top = "10px";
	treeStyleDiv.style.left = "195px";
	treeStyleDiv.style.zIndex = "2000";
	treeStyleDiv.style.display = "none";
	
	treeStyleButton.id = "treeStyleButton";
	treeStyleButton.innerHTML = "Style";
	treeStyleButton.style.width = "80px";
	treeStyleButton.style.height = INSAPHYLOGEO_SIZE_MENU_LABEL_HEIGTH;
	treeStyleButton.style.cursor = "pointer";
	treeStyleButton.style.backgroundColor = "black";
	treeStyleButton.style.borderColor = "white";
	treeStyleButton.style.color = "white";
	treeStyleButton.style.borderRadius = "10px";
	treeStyleButton.style.fontSize = INSAPHYLOGEO_SIZE_FONTSIZE_MENU;
	
	treeStyleDivInside.id = "treeStyleDivInside";
	treeStyleDivInside.class = "dropdown-content";
	treeStyleDivInside.style.display = "none";
	treeStyleDivInside.style.position = "absolute";
	treeStyleDivInside.style.boxShadow = "0 8px 16px 0 rgba(0,0,0,0.9)";
	treeStyleDivInside.style.fontSize = INSAPHYLOGEO_SIZE_FONTSIZE_ITENS_MENU;
	treeStyleDivInside.style.height = "300px";
	treeStyleDivInside.style.overflow = "auto";
	treeStyleDivInside.style.whiteSpace = "nowrap";
	treeStyleDivInside.style.borderRadius = "10px";
	treeStyleDivInside.style.backgroundColor = "black";
	treeStyleDivInside.style.color = "white";
	treeStyleDivInside.style.width = INSAPHYLOGEO_SIZE_MENU_WITH;
	
	propertiesStyleDiv.id = "propertiesStyleDiv";
	propertiesStyleDiv.style.display = "block";
	propertiesStyleDiv.style.paddingLeft = "10px";
	propertiesStyleDiv.style.paddingRight = "10px";
//	propertiesStyleDiv.style.padding = "10px";
	propertiesStyleDiv.style.backgroundColor = "white";
	propertiesStyleDiv.style.color = "black";
	propertiesStyleDiv.style.textAlign = "center";
	
	resetStyleButton.id = "resetStyleButton";
	resetStyleButton.className = "insaphylo_style_buttons_top_menu";
	resetStyleButton.src = insaphylogeo_source_images + "reset.png";
	
	blockStyleButton.id = "blockStyleButton";
	blockStyleButton.className = "insaphylo_style_buttons_top_menu";
	blockStyleButton.src = insaphylogeo_source_images + "non_block.png";
	
	selectNodeColorLabel.id = "selectNodeColorLabel";
	selectNodeColorLabel.innerHTML = "Node Color";
	selectNodeColorLabel.className = "insaphylo_style_item_menu_text";
	
	selectNodeColor.id = "selectNodeColor";
	selectNodeColor.style.display = "block";
	selectNodeColor.style.width = "130px";
	selectNodeColor.style.margin = "auto";
	selectNodeColor.style.cursor = "pointer";
	
	nodeSizeLabel.id = "nodeSizeLabel";
	nodeSizeLabel.innerHTML = "Node Size: 10px";		//  Text
	nodeSizeLabel.className = "insaphylo_style_item_menu_text";
	
	nodeSize.id = "nodeSize";
	nodeSize.type = "range";
	nodeSize.className = "insaphylo_style_item_menu";
	nodeSize.value = "10";
	
	textSizeLabel.id = "textSizeLabel";
	textSizeLabel.innerHTML = "Label Size: 4px";		//  Text
	textSizeLabel.className = "insaphylo_style_item_menu_text";
	
	textSize.id = "textSize";
	textSize.type = "range";
	textSize.value = "12";
	textSize.className = "insaphylo_style_item_menu"
	
	lineWidthLabel.id = "lineWidthLabel";
	lineWidthLabel.innerHTML = "Line Width: 1px";		//  Text
	lineWidthLabel.className = "insaphylo_style_item_menu_text";
	
	lineWidth.id = "lineWidth";
	lineWidth.type = "range";
	lineWidth.value = "20";
	lineWidth.className = "insaphylo_style_item_menu"

	document.getElementById("phylocanvas").appendChild(treeStyleDiv);
	document.getElementById("treeStyleDiv").appendChild(treeStyleButton);
	document.getElementById("treeStyleDiv").appendChild(treeStyleDivInside);
	document.getElementById("treeStyleDivInside").appendChild(propertiesStyleDiv);
	document.getElementById("propertiesStyleDiv").appendChild(resetStyleButton);
	document.getElementById("propertiesStyleDiv").appendChild(blockStyleButton);
	document.getElementById("treeStyleDivInside").appendChild(selectNodeColorLabel);
	document.getElementById("treeStyleDivInside").appendChild(selectNodeColor);
	document.getElementById("treeStyleDivInside").appendChild(nodeSizeLabel);
	document.getElementById("treeStyleDivInside").appendChild(nodeSize);
	document.getElementById("treeStyleDivInside").appendChild(textSizeLabel);
	document.getElementById("treeStyleDivInside").appendChild(textSize);
	document.getElementById("treeStyleDivInside").appendChild(lineWidthLabel);
	document.getElementById("treeStyleDivInside").appendChild(lineWidth);
	
	// -------------------------------------- Metadata --------------------------------------
	
	let treeMetadataDiv = document.createElement("DIV");
	let treeMetadataButton = document.createElement("BUTTON");
	let treeMetadataDivInside = document.createElement("DIV");
	
	let metadataSwitchLabel = document.createElement("LABEL");
	let metadataSwitch = document.createElement("INPUT");
	
	let metadataSwitchLabelDisplay = document.createElement("LABEL");
	let metadataSwitchDisplay = document.createElement("INPUT");
	
	treeMetadataDiv.id = "treeMetadataDiv";
	treeMetadataDiv.classList.add("toggle");
	treeMetadataDiv.style.position = "absolute";
	treeMetadataDiv.style.left = "100px";
	treeMetadataDiv.style.top = "10px";
	treeMetadataDiv.style.zIndex = "2000";
	treeMetadataDiv.style.display = "none";
	
	treeMetadataButton.id = "treeMetadataButton";
	treeMetadataButton.innerHTML = "Metadata";
	treeMetadataButton.style.width = "90px";
	treeMetadataButton.style.height = INSAPHYLOGEO_SIZE_MENU_LABEL_HEIGTH;
	treeMetadataButton.style.cursor = "pointer";
	treeMetadataButton.style.backgroundColor = "black";
	treeMetadataButton.style.borderColor = "white";
	treeMetadataButton.style.color = "white";
	treeMetadataButton.style.borderRadius = "10px";
	treeMetadataButton.style.fontSize = INSAPHYLOGEO_SIZE_FONTSIZE_MENU;
	
	treeMetadataDivInside.id = "treeMetadataDivInside";
	treeMetadataDivInside.class = "dropdown-content";
	treeMetadataDivInside.style.display = "none";
	treeMetadataDivInside.style.position = "absolute";
	treeMetadataDivInside.style.boxShadow = "0 8px 16px 0 rgba(0,0,0,0.9)";
	treeMetadataDivInside.style.fontSize = INSAPHYLOGEO_SIZE_FONTSIZE_ITENS_MENU;
	treeMetadataDivInside.style.height = "250px";
	treeMetadataDivInside.style.overflow = "auto";
	treeMetadataDivInside.style.whiteSpace = "nowrap";
	treeMetadataDivInside.style.borderRadius = "10px";
	treeMetadataDivInside.style.backgroundColor = "black";
	treeMetadataDivInside.style.color = "white";
	treeMetadataDivInside.style.width = INSAPHYLOGEO_SIZE_MENU_WITH;
	
	metadataSwitchLabel.id = "metadataSwitchLabel";
	metadataSwitchLabel.innerHTML = "Select All";
	metadataSwitchLabel.style.display = "block";
	metadataSwitchLabel.style.paddingLeft = "10px";
	metadataSwitchLabel.style.paddingRight = "10px";
	metadataSwitchLabel.style.backgroundColor = "white";
	metadataSwitchLabel.style.color = "black";
	metadataSwitchLabel.style.cursor = "pointer";
	
	metadataSwitch.id = "metadataSwitch";
	metadataSwitch.type = "checkbox";
	metadataSwitch.style.cursor = "pointer";
	metadataSwitch.style.marginLeft = "10px";
	
	metadataSwitchLabelDisplay.id = "metadataSwitchLabelDisplay";
	metadataSwitchLabelDisplay.innerHTML = "Metadata Labels";
	metadataSwitchLabelDisplay.style.display = "block";
	metadataSwitchLabelDisplay.style.paddingLeft = "10px";
	metadataSwitchLabelDisplay.style.paddingRight = "10px";
	metadataSwitchLabelDisplay.style.fontSize = "12px";
	metadataSwitchLabelDisplay.style.backgroundColor = "white";
	metadataSwitchLabelDisplay.style.color = "black";
	metadataSwitchLabelDisplay.style.cursor = "pointer";
	
	metadataSwitchDisplay.id = "metadataSwitchDisplay";
	metadataSwitchDisplay.type = "checkbox";
	metadataSwitchDisplay.style.cursor = "pointer";
	metadataSwitchDisplay.style.marginLeft = "10px";
	
	document.getElementById("phylocanvas").appendChild(treeMetadataDiv);
	
	document.getElementById("treeMetadataDiv").appendChild(treeMetadataButton);
	document.getElementById("treeMetadataDiv").appendChild(treeMetadataDivInside);
	
	document.getElementById("treeMetadataDivInside").appendChild(metadataSwitchLabelDisplay);
	document.getElementById("metadataSwitchLabelDisplay").appendChild(metadataSwitchDisplay);
	
	document.getElementById("treeMetadataDivInside").appendChild(metadataSwitchLabel);
	document.getElementById("metadataSwitchLabel").appendChild(metadataSwitch);
	
	// -------------------------------------- Legend --------------------------------------
	
	let treeLegendDiv = document.createElement("DIV");
	let treeLegendButton = document.createElement("BUTTON");
	let treeLegendDivInside = document.createElement("DIV");
	
	treeLegendDiv.id = "treeLegendDiv";
	treeLegendDiv.classList.add("toggle");
	treeLegendDiv.style.position = "absolute";
	treeLegendDiv.style.left = "10px";
	treeLegendDiv.style.top = "10px";
	treeLegendDiv.style.zIndex = "2000";
	treeLegendDiv.style.display = "none";
	
	treeLegendButton.id = "treeLegendButton";
	treeLegendButton.innerHTML = "Legend";
	treeLegendButton.style.width = "85px";
	treeLegendButton.style.height = INSAPHYLOGEO_SIZE_MENU_LABEL_HEIGTH;
	treeLegendButton.style.cursor = "pointer";
	treeLegendButton.style.backgroundColor = "black";
	treeLegendButton.style.borderColor = "white";
	treeLegendButton.style.color = "white";
	treeLegendButton.style.borderRadius = "10px";
	treeLegendButton.style.fontSize = INSAPHYLOGEO_SIZE_FONTSIZE_MENU
	
	treeLegendDivInside.id = "treeLegendDivInside";
	treeLegendDivInside.class = "dropdown-content";
	treeLegendDivInside.style.display = "none";
	treeLegendDivInside.style.position = "absolute";
	treeLegendDivInside.style.boxShadow = "0 8px 16px 0 rgba(0,0,0,0.9)";
	treeLegendDivInside.style.fontSize = INSAPHYLOGEO_SIZE_FONTSIZE_ITENS_MENU;
	treeLegendDivInside.style.height = "300px";
	treeLegendDivInside.style.overflow = "auto";
	treeLegendDivInside.style.whiteSpace = "nowrap";
	treeLegendDivInside.style.borderRadius = "10px";
	treeLegendDivInside.style.backgroundColor = "black";
	treeLegendDivInside.style.color = "white";
	treeLegendDivInside.style.width = INSAPHYLOGEO_SIZE_MENU_WITH;
	treeLegendDivInside.style.resize = "vertical";
	
	document.getElementById("phylocanvas").appendChild(treeLegendDiv);
	
	document.getElementById("treeLegendDiv").appendChild(treeLegendButton);
	document.getElementById("treeLegendDiv").appendChild(treeLegendDivInside);
	

	// -------------------------------------- Tree Picker --------------------------------------
	let shapes = ["radial", "rectangular", "diagonal", "hierarchical", "circular"];

	for (let i = 0; i < shapes.length; i++) {

	    document.getElementById(shapes[i]).addEventListener("click", function () {

	            document.getElementById("treeButton").src = insaphylogeo_source_images + shapes[i] + ".png";

	            window.PhylocanvasTree.baseNodeSize = "2";

	            window.PhylocanvasTree.setTreeType(shapes[i]); // Choosing type of tree: takes radial, rectangular, circular, diagonal and hierarchy.
	            window.PhylocanvasTree.setTextSize(document.getElementById("textSize").value);
	            window.PhylocanvasTree.setNodeSize(document.getElementById("nodeSize").value);

	            treeType = shapes[i];

	    });
	}

	// -------------------------------------- Mouse Input Dynamics --------------------------------------

	// -------------------------------- Node Size

	document.getElementById("nodeSize").addEventListener("input", function () {
		window.PhylocanvasTree.baseNodeSize = "2";
		window.PhylocanvasTree.setNodeSize(document.getElementById("nodeSize").value);
	    document.getElementById("nodeSizeLabel").innerHTML = "Node Size: " + document.getElementById("nodeSize").value + "px";
	    document.getElementById("resetStyleButton").style.transform = "rotate(0deg)";
	});

	// -------------------------------- Label Size

	document.getElementById("textSize").addEventListener("input", function () {

	    //window.PhylocanvasTree.draw();
	    //window.PhylocanvasTree.setTreeType(treeType);
		window.PhylocanvasTree.setTextSize(document.getElementById("textSize").value);
	    document.getElementById("textSizeLabel").innerHTML = "Label Size: " + document.getElementById("textSize").value + "px";
	    document.getElementById("resetStyleButton").style.transform = "rotate(0deg)";
	});

	// -------------------------------- Line Width

	document.getElementById("lineWidth").addEventListener("input", function () {

		window.PhylocanvasTree.lineWidth = document.getElementById("lineWidth").value/20;

		window.PhylocanvasTree.setNodeSize(document.getElementById("nodeSize").value);
		window.PhylocanvasTree.setTextSize(document.getElementById("textSize").value);

	    document.getElementById("lineWidthLabel").innerHTML = "Line Width: " + document.getElementById("lineWidth").value/20 + "px";
	    document.getElementById("resetStyleButton").style.transform = "rotate(0deg)";
	});

	// -------------------------------------- Mouse Over Dynamics --------------------------------------

	// -------------------------------- Mouse Over

	document.getElementById("treeDivInside").addEventListener("mouseover", function () {

	    treeDivInside.style.display = "block";
	});

	document.getElementById("treeButton").addEventListener("mouseover", function () {

	    treeDivInside.style.display = "block";
	});

	document.getElementById("treeStyleDivInside").addEventListener("mouseover", function () {

	    treeStyleDivInside.style.display = "block";
	});

	document.getElementById("treeStyleButton").addEventListener("mouseover", function () {

	    treeStyleDivInside.style.display = "block";
	});

	document.getElementById("treeMetadataDivInside").addEventListener("mouseover", function () {

	    treeMetadataDivInside.style.display = "block";
	});

	document.getElementById("treeMetadataButton").addEventListener("mouseover", function () {

	    treeMetadataDivInside.style.display = "block";
	});

	document.getElementById("treeLegendDivInside").addEventListener("mouseover", function () {

	    treeLegendDivInside.style.display = "block";
	});

	document.getElementById("treeLegendButton").addEventListener("mouseover", function () {

	    treeLegendDivInside.style.display = "block";
	});

	// -------------------------------- Mouse Out

	document.getElementById("treeDivInside").addEventListener("mouseout", function () {

	    treeDivInside.style.display = "none";
	});

	document.getElementById("treeButton").addEventListener("mouseout", function () {

	    treeDivInside.style.display = "none";
	});

	document.getElementById("treeStyleDivInside").addEventListener("mouseout", function () {

	    treeStyleDivInside.style.display = "none";
	});

	document.getElementById("treeStyleButton").addEventListener("mouseout", function () {

	    if (document.getElementById("blockStyleButton").src.search(insaphylogeo_source_images + "non_block.png") !== -1) {

	        treeStyleDivInside.style.display = "none";
	    }
	});

	document.getElementById("treeMetadataDivInside").addEventListener("mouseout", function () {

	    treeMetadataDivInside.style.display = "none";
	});

	document.getElementById("treeMetadataButton").addEventListener("mouseout", function () {

	    treeMetadataDivInside.style.display = "none";
	});

	document.getElementById("treeLegendDivInside").addEventListener("mouseout", function () {

	        treeLegendDivInside.style.display = "none";
	});

	document.getElementById("treeLegendButton").addEventListener("mouseout", function () {

		var element = document.getElementById("block")
	    if (element == null){
	    	treeLegendDivInside.style.display = "none";
	    }
	    else if (element.src.search(insaphylogeo_source_images + "non_block.png") !== -1) {

	        treeLegendDivInside.style.display = "none";
	    }
	});
}
// -------------------------------------- Toggle Wheel --------------------------------------

function toggle() {

    let x = document.getElementsByClassName("toggle");

    for (let i = 0; i < x.length; i++) {
        if (x[i].style.display === "none") {
            x[i].style.display = "block";
            document.getElementById("treeButton").style.display = "block";
            document.getElementById("toggle").style.transform = "rotate(-15deg)";
        } else {
            x[i].style.display = "none";
            document.getElementById("treeButton").style.display = "none";
            document.getElementById("toggle").style.transform = "rotate(15deg)";
        }
    }
}


//-------------------------------------- Select All Button --------------------------------------
function selectAlll (metaData, max) {

    document.getElementById("metadataSwitch").addEventListener("change", function () {

        if (document.getElementById("metadataSwitch").checked === true) {

            for (let i = 0; i < document.getElementsByClassName("metadataSwitch").length; i++) {

                document.getElementsByClassName("metadataSwitch")[i].checked = true;

            }

            for (let j = 0; j < max; j++) {

                        for (let i = 0; i < window.PhylocanvasTree.leaves.length; i++) {

                            (window.PhylocanvasTree.leaves[i].data)[Object.keys(metaData[window.PhylocanvasTree.leaves[i].label])[j]] = (window.PhylocanvasTree.leaves[i].data)[Object.keys(metaData[window.PhylocanvasTree.leaves[i].label])[j]] || {};

                            (window.PhylocanvasTree.leaves[i].data)[Object.keys(metaData[window.PhylocanvasTree.leaves[i].label])[j]]["colour"] = color[Object.keys(metaData[window.PhylocanvasTree.leaves[i].label])[j] + Object.values(metaData[window.PhylocanvasTree.leaves[i].label])[j]];

                            if (document.getElementById("metadataSwitchDisplay").checked === true) {

                                (window.PhylocanvasTree.leaves[i].data)[Object.keys(metaData[window.PhylocanvasTree.leaves[i].label])[j]]["label"] = Object.values(metaData[window.PhylocanvasTree.leaves[i].label])[j] || "No Data";

                            } else {

                                (window.PhylocanvasTree.leaves[i].data)[Object.keys(metaData[window.PhylocanvasTree.leaves[i].label])[j]]["label"] = "";

                            }
                        }

                window.PhylocanvasTree.setTreeType(treeType); // Choosing type of tree: takes radial, rectangular, circular, diagonal and hierarchy.
                window.PhylocanvasTree.setNodeSize(document.getElementById("nodeSize").value);
                window.PhylocanvasTree.setTextSize(document.getElementById("textSize").value);
            }

        } else {

            for (let i = 0; i < document.getElementsByClassName("metadataSwitch").length; i++) {

                document.getElementsByClassName("metadataSwitch")[i].checked = false;

            }

                for (let i = 0; i < window.PhylocanvasTree.leaves.length; i++) {

                    window.PhylocanvasTree.leaves[i].data = {
                    };

                }

            window.PhylocanvasTree.setTreeType(treeType); // Choosing type of tree: takes radial, rectangular, circular, diagonal and hierarchy.
            window.PhylocanvasTree.setNodeSize(document.getElementById("nodeSize").value);
            window.PhylocanvasTree.setTextSize(document.getElementById("textSize").value);

        }
    });
}

// -------------------------------------- Metadata Labels Button --------------------------------------

function displayLabel (metaData, max) {

    document.getElementById("metadataSwitchDisplay").addEventListener("change", function () {

    let metaLabels = [];

    for (let i = 0; i < document.getElementsByClassName("metadataSwitch").length; i++) {

        if (document.getElementsByClassName("metadataSwitch")[i].checked === true) {

            metaLabels.push(document.getElementsByClassName("metadataSwitchLabel")[i].innerText);

        }
    }

    for (let j = 0; j < max; j++) {

                for (let i = 0; i < window.PhylocanvasTree.leaves.length; i++) {

                    if (metaLabels.find(function(element) {
                        return (element === Object.keys(metaData[window.PhylocanvasTree.leaves[i].label])[j]);
                    }) !== undefined) {

                        (window.PhylocanvasTree.leaves[i].data)[Object.keys(metaData[window.PhylocanvasTree.leaves[i].label])[j]] = (window.PhylocanvasTree.leaves[i].data)[Object.keys(metaData[window.PhylocanvasTree.leaves[i].label])[j]] || {};

                        (window.PhylocanvasTree.leaves[i].data)[Object.keys(metaData[window.PhylocanvasTree.leaves[i].label])[j]]["colour"] = color[Object.keys(metaData[window.PhylocanvasTree.leaves[i].label])[j] + Object.values(metaData[window.PhylocanvasTree.leaves[i].label])[j]];

                        if (document.getElementById("metadataSwitchDisplay").checked === true) {

                            (window.PhylocanvasTree.leaves[i].data)[Object.keys(metaData[window.PhylocanvasTree.leaves[i].label])[j]]["label"] = Object.values(metaData[window.PhylocanvasTree.leaves[i].label])[j] || "No Data";

                        } else {

                            (window.PhylocanvasTree.leaves[i].data)[Object.keys(metaData[window.PhylocanvasTree.leaves[i].label])[j]]["label"] = "";

                        }
                    }
                }

        window.PhylocanvasTree.setTreeType(treeType);
        window.PhylocanvasTree.setNodeSize(document.getElementById("nodeSize").value);
        window.PhylocanvasTree.setTextSize(document.getElementById("textSize").value);
    }
    })
}