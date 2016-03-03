#!/usr/bin/env node
var argv = require('minimist')(process.argv.slice(2));

var url = argv.swagger || 'http://data.gramene.org/v50/swagger';
var feature = argv.feature || 'gene';
var combiner = argv.combiner || 'canonical';
var idFile = argv.ids;
var batchSize = argv.batchSize || 1000;

global.gramene = {defaultServer: url};
var gramene = require('gramene-search-client').client.grameneClient;
if (idFile) {
  var ids = [];
  require('readline').createInterface({
    input: require('fs').createReadStream(idFile),
    terminal: false
  })
  .on('line', function(line) {
    ids.push(line);
  })
  .on('close', function() {
    gramene.then(function(client) {
      for(var i=0; i<ids.length; i+=batchSize) {
        var params = {
          rows: batchSize,
          wt: 'bed',
          bedFeature: feature,
          bedCombiner: combiner,
          idList: ids.slice(i, i+batchSize)
        };
        if (argv.taxon_id) {
          params.taxon_id = argv.taxon_id;
        }
        client['Data access'].genes(params, function(res) {
          console.log(res.data);
        })
      }
    })
  })
}
else if (argv.taxon_id) {
  gramene.then(function(client) {
    client['Data access'].genes(
      {
        wt: 'bed',
        rows: -1,
        bedFeature: feature,
        bedCombiner: combiner,
        taxon_id: argv.taxon_id,
        db_type: 'core'
      }, function(res) {
        console.log(res.data);
      }
    )
  })
}
