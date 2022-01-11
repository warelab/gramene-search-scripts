#!/usr/bin/env node
var through2 = require('through2');
var byline = require('byline');
var fs = require('fs');
var argv = require('minimist')(process.argv.slice(2));

var url = argv.swagger || 'https://data.gramene.org/vitis1/swagger';
var idFile = argv.ids;
var vitis_taxa = [
	297600000,
	297600001,
	297600002,
	297600003,
	297600004,
	297600005,
	297600006,
	297600007,
	297600008,
	297600009
];
var sbi_taxa = [
	4558,
	1000655996,
	1000561071,
	1000656001,
	1000651496
];
var taxa = sbi_taxa;
global.gramene = {defaultServer: url};
var gramene = require('gramene-search-client').client.grameneClient;
let GrameneTrees = require('gramene-trees-client');

var reader = byline(fs.createReadStream(idFile));

var treeFetcher = through2.obj(function(id, enc, done) {
  var that = this;
  id = id.toString();
  gramene.then(function(client) {
    client['Data access']['genetrees']({idList:[id],rows:-1}).then(function(res) {
      that.push(res.obj);
      done();
    })
  })
});

function walkTo(node, test, process) {
  if (test(node)) {
    process(node);
  }
  else {
    if (node.hasChildren()) {
      node.children.forEach(function(childNode) {
        walkTo(childNode, test, process)
      })
    }
  }
}

var checkTree = function checkTree() {
  var results = [];
	var cols = ['tree','taxa'].concat(taxa);
	results.push(cols.join("\t"));
  var transform = function(tree, enc, done) {
    let genetree = GrameneTrees.genetree.tree(tree);
		let tally = {};
		taxa.forEach(tid => tally[tid]=0);
		walkTo(
			genetree,
			node => (!node.children || node.children.length === 0) && tally.hasOwnProperty(node.model.taxon_id),
			node => tally[node.model.taxon_id]++
		);
		let row = [];
		row.push(genetree._id);
		row.push(0); // number of taxa present
		for (const tid in taxa) {
			row.push(tally[tid]);
			if (tally[tid] > 0) {
				row[1]++;
			}
		}
		if (row[1] > 0) {
			results.push(row.join("\t"));
		}
    done();
  };
  
  var flush = function(done) {
	  this.push(results.join("\n"));
    done();
  }
  
  return through2.obj(transform, flush);
};


reader
.pipe(treeFetcher)
.pipe(checkTree())
.pipe(process.stdout);
