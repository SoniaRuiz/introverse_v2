
function loadBigWigFiles(bigwigCategories, bigwigURLs, bigwigTypes) {
  
  // Destroy any previous igv instances
  igv.removeAllBrowsers();
  
  // GET THE GENOMIC COORDINATES TO VISUALISE
  var igvDiv = document.getElementById("igv-div");
  var chr = document.getElementById("chr_tab1").value;
  var start = document.getElementById("start_tab1").value;
  var end = document.getElementById("end_tab1").value;
  var strand = document.getElementById("strand_tab1").value;
  
  var genomic_locus = "pten";
  if (chr !== "") {
    genomic_locus = "chr" + chr + ":" + start + "-" + end;
  }
  
  var reference = {
    tracks:[]
  };
  reference.tracks.push({
    "name": "Refseq Genes",
    "url": "https://s3.amazonaws.com/igv.org.genomes/hg38/refGene.txt.gz",
    "order": 1000000,
    "height": 200,
    "indexed": false
    });
          


  // GET THE BIGWIG categories to process
  var categories = bigwigCategories.split(';');
  let uniqueCategories = [...new Set(categories)];
  
  
  // GET THE BIGWIG types TO LOAD
  var strand_type = bigwigTypes.split(';');
  
  // GET THE BIGWIG URLs TO LOAD
  var URLs = bigwigURLs.split(';');
  

  for(var i = 0; i < uniqueCategories.length; i++)
  {
    
    var max_tracks = "";
    
    if (URLs.length < 50) {
      max_tracks = URLs.length;
    } else if (uniqueCategories.length == 1 && URLs.length > 50) {
      max_tracks = 20;
    } else if (uniqueCategories.length > 1 && ((URLs.length/uniqueCategories.length) > 100)) {
      max_tracks = 20 * uniqueCategories.length;
    } else {
      max_tracks = URLs.length;
    }
    
    var track_name = uniqueCategories[i];
    if ( uniqueCategories[i] == "case" ) {
      track_name = "shRNA knockdown";
    } 
    var local_bigwig_track = {
      name: track_name,
      type: "merged",
      height: 100,
      tracks: []
    };
    
    
     
    for(var j = 0; j < max_tracks; j++)
    {
      
      if (uniqueCategories[i] == categories[j]) {
        
        var coverage_color = "#999999";
        
        if ( strand_type[j] == "minus" ) {
          coverage_color = "#cccccc";
        }
        
        local_bigwig_track.tracks.push({ 
            "type" : 'wig',
            "format" : "bigwig",
            "url" : URLs[j], 
            "color": coverage_color
        });
      }
      
    }
    
    
    local_bigwig_track.name = local_bigwig_track.name + " (displaying " + local_bigwig_track.tracks.length + " overlaid BigWig samples)";
    reference.tracks.push(local_bigwig_track);
    console.log(URLs.length);
  }
  
  
  console.log( JSON.stringify(reference) );

  
  // CONFIGURE THE IGV BROWSER
  var options =
  {
    reference: {
      "id": "hg38",
      "name": "Human (GRCh38/hg38)",
      
      "fastaURL": "https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg38/hg38.fa",
      "indexURL": "https://s3.amazonaws.com/igv.broadinstitute.org/genomes/seq/hg38/hg38.fa.fai",
      "cytobandURL": "https://s3.amazonaws.com/igv.broadinstitute.org/annotations/hg38/cytoBandIdeo.txt",
      "tracks": reference.tracks
    },
    locus: genomic_locus
  };

  console.log( JSON.stringify(options) );
  
  igv.createBrowser(igvDiv, options)
          .then(function (browser) {
              console.log("Created IGV browser!");
          });
  
  
}

