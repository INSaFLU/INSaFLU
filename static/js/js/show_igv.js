//
//	This test on the name change
//
//$("#phylocanvas").on("change paste keyup", function () {

$('#collapseOne').on('shown.bs.collapse', function () {
    show_igv();
});

//remove the tree from the screen
$('#collapseOne').on('hidden.bs.collapse', function () {
    $('#show_igv_id').empty();
    $('#bam_file_id').empty();
});

// https://github.com/igvteam/igv.js/issues/241
// https://developer.mozilla.org/en-US/docs/Web/HTTP/CORS
// https://github.com/racitup/static-ranges

// show igv
function show_igv() {
    
    $.ajax({
        /// spin 
        beforeSend: function() {
            $('#loader_igv_id').show();
        },
        complete: function(){
            $('#loader_igv_id').hide();
        },
        
        data : { 
            project_sample_id : $('#show_igv_id').attr("project_sample_id"),
        }, // data sent with the get request
        
        url: $('#show_igv_id').attr("show-igv-url"),
        success: function (data) {
            if (data['is_ok']) {
                
                /// set the files names
                $('#bam_file_id').empty()
                $('#bam_file_id').append(data['bam_file_id'])
                $('#bai_file_id').empty()
                $('#bai_file_id').append(data['bai_file_id'])
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
                browser = igv.createBrowser(document.getElementById("show_igv_id"), options);
            }
            else{
                $('#show_igv_id').append('<div class="alert alert-warning alert-dismissable"><strong>Fail</strong> to load bam file.</div>')
            }
        },
        
        // handle a non-successful response
        error : function(xhr,errmsg,err) {
            $('#show_igv_id').append('<div class="alert alert-warning alert-dismissable"><strong>Fail</strong> to load bam file.</div>')
            console.log(xhr.status + ": " + xhr.responseText); // provide a bit more info about the error to the console
        }
    });
}