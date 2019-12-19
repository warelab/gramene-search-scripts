#!/usr/bin/env node
var through2 = require('through2');
var byline = require('byline');
var fs = require('fs');
var argv = require('minimist')(process.argv.slice(2));

var url = argv.swagger || 'http://data.gramene.org/maize3/swagger';
var idFile = argv.ids;

global.gramene = {defaultServer: url};
var gramene = require('gramene-search-client').client.grameneClient;
let GrameneTrees = require('gramene-trees-client');

var reader = byline(fs.createReadStream(idFile));

var fetcher = through2.obj(function(id, enc, done) {
  var that = this;
  id = id.toString();
  gramene.then(function(client) {
    client['Data access']['genetrees']({idList:[id],rows:-1}).then(function(res) {
      that.push(res.obj);
      done();
    })
  })
});


var checkTree = function checkTree() {
  var results = [];
  var transform = function(tree, enc, done) {
    let genetree = GrameneTrees.genetree.tree(tree);
    GrameneTrees.extensions.addConsensus(genetree);
    // find Zea (taxon_id 4575) speciation nodes (bfs traversal)
    // compare coverage of each protein in the subtree to the consensus at this node
    // for each protein output gene id, protein id, taxon_id, and the ratio above and pct identity with consensus
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
    
    function isZeaSpeciation(node) {
      return (node.model.taxon_name === "Zea" && node.model.node_type === "speciation")
    }
    function isAndropogoneaeSpeciation(node) {
      return (node.model.taxon_name === "Andropogoneae" && node.model.node_type === "speciation")
    }
    function isPoaceaeSpeciation(node) {
      return (node.model.taxon_name === "Poaceae" && node.model.node_type === "speciation")
    }
    
    function coverage_similarity(a,b) {
      let seqA = a.model.consensus.sequence;
      let seqB = b.model.consensus.sequence;
      let totalA=0;
      let totalB=0;
      let aligned=0;
      const gapCode = '-'.charCodeAt(0);
      for(var i=0; i<seqA.length; i++) {
        if (seqA[i] !== gapCode) {
          totalA++;
          if (seqB[i] !== gapCode) {
            aligned++;
          }
        }
        if (seqB[i] !== gapCode) {
          totalB++;
        }
      }
      return (aligned/totalA + aligned/totalB)/2;
    }

    function compareToConsensus(node) {
      node.all(function(leaf) {
        if (!leaf.hasChildren()) {
          
          function coverage(a,b) {
            let seqA = a.model.consensus.sequence;
            let seqB = b.model.consensus.sequence;
            let total=0;
            let aligned=0;
            const gapCode = '-'.charCodeAt(0);
            for(var i=0; i<seqA.length; i++) {
              if (seqA[i] !== gapCode) {
                total++;
                if (seqB[i] !== gapCode) {
                  aligned++;
                }
              }
            }
            return aligned/total;
          }
          results.push([
            genetree.model.tree_stable_id,
            node.model.node_id,
            node.model.consensus.nSeqs,
            leaf.model.taxon_name,
            coverage(node, leaf),
            GrameneTrees.extensions.identity(leaf, node),
            leaf.model.gene_stable_id,
            leaf.model.protein_stable_id
          ].join("\t"))
        }
      })
    }

    function compareToConsensusOfSorghum(nodeA) {
      let sorghum_node;
      let maize_node;
      nodeA.children.forEach(function(childNode) {
        if (childNode.model.taxon_name === "Sorghum bicolor") {
          sorghum_node = childNode;
        }
        else {
          maize_node = childNode;
        }
      });
      if (maize_node && sorghum_node) {
        maize_node.all(function(leaf) {
          if (!leaf.hasChildren()) {
            results.push([
              genetree.model.tree_stable_id,
              nodeA.model.node_id,
              nodeA.model.consensus.nSeqs,
              leaf.model.taxon_name,
              coverage_similartiy(sorghum_node, leaf),
              GrameneTrees.extensions.identity(leaf, sorghum_node),
              leaf.model.gene_stable_id
            ].join("\t"))
          }
        })
      }
    }
    
    function compareSorghumToConsensusOfRice(poaceae_node) {
      let rice_node;
      walkTo(poaceae_node,
        function(node) {
          return (node.model.taxon_id === 39947)
        },
        function(node) {
          rice_node = node;
        }
      );
      let sorghum_node;
      walkTo(poaceae_node,
        function(node) {
          return (node.model.taxon_id === 4558)
        },
        function(node) {
          sorghum_node = node;
        }
      );
      
      if (rice_node && sorghum_node) {
        sorghum_node.all(function(leaf) {
          if (!leaf.hasChildren()) {
            results.push([
              genetree.model.tree_stable_id,
              poaceae_node.model.node_id,
              poaceae_node.model.consensus.nSeqs,
              leaf.model.taxon_name,
              coverage_similarity(rice_node, leaf),
              leaf.model.gene_stable_id
            ].join("\t"))
          }
        })
      }
    }

    walkTo(genetree, isPoaceaeSpeciation, compareSorghumToConsensusOfRice);
    done();
  };
  
  var flush = function(done) {
    this.push(results.join("\n"));
    done();
  }
  
  return through2.obj(transform, flush);
};


reader
.pipe(fetcher)
.pipe(checkTree())
.pipe(process.stdout);
