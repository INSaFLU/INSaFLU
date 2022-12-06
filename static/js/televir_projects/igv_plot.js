
/* televir_projects/result_display.html
* This script is used to display an IGV plot of a bam file
*/

$(document).on('click', '.btn', function (e) {
console.log("igv_browse");
if (e.target.id == "igv_browse") 
    var accid=$(this).attr('accid');
    console.log('igv_display_'+accid);
    var igv_display = document.getElementById('igv_display_' + accid);
    var igv_display_status= igv_display.style.display;
    var igv_display_className= igv_display.className;
    console.log(igv_display_status);
    console.log(igv_display_className);
    console.log(accid);
    console.log(/\bopen\b/.test(igv_display_className)); 
    console.log(igv_display.className.replace(" open",''));   
    
    if (/\bopen\b/.test(igv_display_className)){
        igv_display.className= igv_display.className.replace(" open",'');
        var show_igv_div = document.getElementById('show_igv_' + accid);
        document.getElementById('show_igv_' + accid).innerHTML = ""

    }else{
        igv_display.className += " open";
        show_igv($(this));
    }
});



function show_igv(item) {
        var accid=item.attr('accid');
        var project_pk= item.attr('project_pk');
        var sample_pk= item.attr('sample_pk');
        var run_pk = item.attr('run_pk');
        var unique_id = item.attr('reference_id');
        var url = "{% url 'igv_browser' %}";
        console.log("show_igv");
        console.log(url);

        $.ajax({
            /// spin 
            beforeSend: function() {
                $('#igv_display_' + accid).show();
            },
            complete: function(){
                $('#igv_display_' + accid).hide();
            },
            
            data : { 
                'project_pk': project_pk,
                'sample_pk': sample_pk,
                'run_pk': run_pk,
                'unique_id': unique_id,
                'accid': accid,
            }, // data sent with the get request
            
            url: url,

            success: function (data) {
                if (data['is_ok']) {
                    
                    /// set the files names
                    $('#bam_file_id').empty()
                    $('#bam_file_id').append(data['bam_file_id'])
                    $('#bai_file_id').empty()
                    $('#bai_file_id').append(data['bai_file_id'])
                    $('#vcf_file_id').empty()
                    $('#vcf_file_id').append(data['vcf_file_id'])
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
                    console.log("options");
                    console.log(document.getElementById('show_igv_' + accid));
                    browser = igv.createBrowser(document.getElementById('show_igv_' + accid), options);
                }
                else{
                    $('#show_igv_' + accid).append('<div class="alert alert-warning alert-dismissable"><strong>Fail</strong> to load bam file.</div>')
                }
            },
            
            // handle a non-successful response
            error : function(xhr,errmsg,err) {
                $('#show_igv_' + accid).append('<div class="alert alert-warning alert-dismissable"><strong>Fail</strong> to load bam file.</div>')
                console.log(xhr.status + ": " + xhr.responseText); // provide a bit more info about the error to the console
            }
        });
    }