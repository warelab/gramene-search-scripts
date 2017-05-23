#!/usr/bin/env node
var through2 = require('through2');
var byline = require('byline');
var fs = require('fs');
var argv = require('minimist')(process.argv.slice(2));

var url = argv.swagger || 'http://data.gramene.org/swagger';
var idFile = argv.ids;
var batchSize = argv.batchSize || 100;

global.gramene = {defaultServer: url};
var gramene = require('gramene-search-client').client.grameneClient;
var ggp = require('gramene-gene-positions');
var _ = require('lodash');

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
      wt: 'json',
      idList: ids
    };
    client['Data access'].genes(params, function(res) {
      var genes = JSON.parse(res.data);
      genes.forEach(function(gene) {
        that.push(gene);
      });
      done();
    })
  })
});

var domainExtractor = through2.obj(function(gene, enc, done) {
  var that = this;
  if (_.has(gene, 'gene_structure.transcripts')) {
    gene.gene_structure.transcripts.forEach(function (transcript) {
      if (_.has(transcript, 'translation.features.domain.entries')) {
        transcript.translation.features.domain.entries.forEach(function (domain) {
          var genomicStart = ggp.remap(gene, domain.start, 'protein', 'genome');
          var genomicEnd = ggp.remap(gene, domain.end, 'protein', 'genome');
          var projectedDomain = _.clone(domain);
          if (genomicStart > genomicEnd) {
            projectedDomain.gstart = genomicEnd;
            projectedDomain.gend = genomicStart;
          }
          else {
            projectedDomain.gstart = genomicStart;
            projectedDomain.gend = genomicEnd;
          }
          projectedDomain.geneId = gene._id;
          projectedDomain.translationId = transcript.translation.id;
          projectedDomain.region = gene.location.region;
          that.push(projectedDomain);
        });
      }
    });
  }
  done();
});

var keyList = ['region','gstart','gend','gene','geneId','translationId','interpro','start','end','db','name','description'];
var tabber = through2.obj(function(domain, enc, done) {
  var vals = keyList.map(function(k) {
    return domain[k];
  });
  this.push(vals.join("\t") + "\n");
  done();
});

console.log(keyList.join("\t"));
reader
  .pipe(batcher(batchSize))
  .pipe(fetcher)
  .pipe(domainExtractor)
  .pipe(tabber)
  .pipe(process.stdout);
