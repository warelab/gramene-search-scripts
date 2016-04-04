#!/usr/bin/env node
var through2 = require('through2');
var byline = require('byline');
var fs = require('fs');
var argv = require('minimist')(process.argv.slice(2));

var url = argv.swagger || 'http://data.gramene.org/v50/swagger';
var idFile = argv.ids;
var facet = argv.facet;
var batchSize = argv.batchSize || 100;

global.gramene = {defaultServer: url};
var gramene = require('gramene-search-client').client.grameneClient;

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
      rows: 0,
      q: 'id:(' + ids.join(' OR ') + ')',
      'facet.field' : '{!facet.limit=-1 facet.mincount=1}' + facet,
      facet: true
    };
    client['Search'].genes(params, function(res) {
      that.push(res.data);
      done();
    })
  })
});

var expanded = {};
var extract = function extract() {

  var transform = function(doc, enc, done) {
    var a = JSON.parse(doc).facet_counts.facet_fields[facet];
    for(var i=0;i<a.length; i+=2) {
      expanded[a[i]] = 1;
    }
    done();
  };
  
  var flush = function(done) {
    this.push(Object.keys(expanded).join("\n"));
    done();
  }
  
  return through2.obj(transform, flush);
};


reader
.pipe(batcher(batchSize))
.pipe(fetcher)
.pipe(extract())
.pipe(process.stdout);
