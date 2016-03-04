#!/usr/bin/env node
var through2 = require('through2');
var byline = require('byline');
var fs = require('fs');
var argv = require('minimist')(process.argv.slice(2));

var url = argv.swagger || 'http://data.gramene.org/v50/swagger';
var feature = argv.feature || 'gene';
var combiner = argv.combiner || 'canonical';
var idFile = argv.ids;
var batchSize = argv.batchSize || 100;

global.gramene = {defaultServer: url};
var gramene = require('gramene-search-client').client.grameneClient;
if (idFile) {
  var reader = byline(fs.createReadStream(idFile));
  var ids = [];

  var batcher = function batcher(batchSize) {
    
    var transform = function (line, enc, done) {
      if (ids.length === batchSize) {
        var batch = ids;
        this.push(batch);
        ids = [];
      }
      ids.push(line);
      done();
    };
    
    var flush = function(done) {
      this.push(ids);
      done();
    }
    
    return through2.obj(transform,flush);
  };

  var fetcher = through2.obj(function(ids, enc, done) {
    var that = this;
    gramene.then(function(client) {
      var params = {
        rows: batchSize,
        wt: 'bed',
        bedFeature: feature,
        bedCombiner: combiner,
        idList: ids
      };
      if (argv.taxon_id) {
        params.taxon_id = argv.taxon_id;
      }
      client['Data access'].genes(params, function(res) {
        that.push(res.data);
        done();
      })
    })
  });

  reader
  .pipe(batcher(batchSize))
  .pipe(fetcher)
  .pipe(process.stdout);
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
