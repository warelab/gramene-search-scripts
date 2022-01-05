#!/usr/bin/env node
var through2 = require('through2');
var byline = require('byline');
var fs = require('fs');
var argv = require('minimist')(process.argv.slice(2));

var url = argv.swagger || 'https://data.gramene.org/vitis1/swagger';
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
    function isMaizeSubtree(node) {
      const isMaizeRE = new RegExp(/zea/i);
      if (isMaizeRE.test(node.model.taxon_name))
        return true;
      if (node.model.taxon_name === "unknown") {
        let isMaize = true;
        node.children.forEach(function(childNode) {
          isMaize &= isMaizeSubtree(childNode)
        });
        return isMaize;
      }
      return false;
    }
    function isSorghumSubtree(node) {
      const isSorghumRE = new RegExp(/sorghum/i);
      if (isSorghumRE.test(node.model.taxon_name))
        return true;
      if (node.model.taxon_name === "unknown") {
        let isSorghum = true;
        node.children.forEach(function(childNode) {
          isSorghum &= isSorghumSubtree(childNode)
        });
        return isSorghum;
      }
      return false;
    }
    function isOryzaSubtree(node) {
      const isOryzaRE = new RegExp(/oryza/i);
      if (isOryzaRE.test(node.model.taxon_name))
        return true;
      if (node.model.taxon_name === "unknown") {
        let isOryza = true;
        node.children.forEach(function(childNode) {
          isOryza &= isOryzaSubtree(childNode)
        });
        return isOryza;
      }
      return false;
    }
    function isNotSorghumSpeciation(node) {
      const isSorghum = new RegExp(/sorghum/i);
      return (node.model.node_type === "speciation" && !isSorghum.test(node.model.taxon_name))
    }
    function isSorghumSpeciation(node) {
      const isSorghum = new RegExp(/sorghum/i);
      return (node.model.node_type === "speciation" && isSorghum.test(node.model.taxon_name))
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
	function isPN(node) {
		return (node.model.taxon_id === 29760 && !node.hasChildren())
	}
	function isFlagged(node) {
		return (!node.hasChildren() && (node.model.taxon_id === 29760 || node.model.taxon_id === 297600000 || node.model.taxon_id === 4558))
	}
	function isRosidSpeciation(node) {
		return (node.model.taxon_name === "rosids" && node.model.node_type === "speciation")
	}
    
    function coverage_similarity(a,b) {
      let seqA = a.model.consensus.sequence;
      let seqB = b.model.consensus.sequence;

      // make a vector of weights that normalize the frequency in the msa
      let weight = [];
      a.model.consensus.frequency.forEach(function(freq) {
         weight.push(freq/a.model.consensus.nSeqs);
      });
      let totalA=0;
      let aligned=0;
      let matches=0;
      const gapCode = '-'.charCodeAt(0);
      for(var i=0; i<seqA.length; i++) {
        if (seqA[i] !== gapCode) {
          totalA += weight[i];

          if (seqB[i] !== gapCode) {
            aligned += weight[i];
            if (seqB[i] === seqA[i]) {
              matches += weight[i];
            }
            else {
              matches -= weight[i];
            }
          }
        }
      }
      return matches <  0 ? 0 : matches/totalA;
    }

    // function coverage_similarity(a,b) {
    //   let seqA = a.model.consensus.sequence;
    //   let seqB = b.model.consensus.sequence;
    //   let total=0;
    //   let aligned=0;
    //   const gapCode = '-'.charCodeAt(0);
    //   for(var i=0; i<seqA.length; i++) {
    //     if (seqA[i] !== gapCode) {
    //       total++;
    //       if (seqB[i] !== gapCode) {
    //         aligned++;
    //       }
    //     }
    //   }
    //   return aligned/total;
    // }

    function compareToConsensus(node) {
      node.all(function(leaf) {
        if (!leaf.hasChildren()) {
          results.push([
            genetree.model.tree_stable_id,
            genetree.model.taxon_id,
            node.model.node_id,
            node.model.consensus.nSeqs,
            leaf.model.taxon_id,
            leaf.model.system_name,
            coverage_similarity(node, leaf),
            leaf.parent.model.node_type,
            leaf.model.gene_stable_id
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

    function compareToConsensusOfArabidopsis(nodeA) {
      let ath_node;
      let grape_node;
      nodeA.children.forEach(function(childNode) {
		  console.log(childNode.model.taxon_id);
        if (childNode.model.taxon_id === 3702) {
          ath_node = childNode;
        }
        else {
          grape_node = childNode;
        }
      });
      if (ath_node && grape_node) {
        grape_node.all(function(leaf) {
          if (!leaf.hasChildren()) {
            results.push([
              genetree.model.tree_stable_id,
              nodeA.model.node_id,
              nodeA.model.consensus.nSeqs,
              leaf.model.taxon_id,
              coverage_similarty(ath_node, leaf),
              GrameneTrees.extensions.identity(leaf, ath_node),
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

	function getSince(grape_gene) {
		let node = grape_gene;
		while (node.parent && node.model.taxon_id !== 3398) {
			node = node.parent;
		}
		results.push([
			grape_gene.model.gene_stable_id,
			node.model.node_id,
			node.model.taxon_id,
			node.model.taxon_name
		].join("\t"));
	}
    // walkTo(genetree, isPoaceaeSpeciation, compareSorghumToConsensusOfRice);
    // walkTo(genetree, isRosidSpeciation, compareToConsensus);
    // walkTo(genetree, isSorghumSpeciation, compareToConsensus);
   // walkTo(genetree, isSorghumSubtree, compareToConsensus);
    // walkTo(genetree, isMaizeSubtree, compareToConsensus);
    walkTo(genetree, isOryzaSubtree, compareToConsensus);
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
