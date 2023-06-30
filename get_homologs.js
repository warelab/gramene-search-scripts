#!/usr/bin/env node
const fetch = require('node-fetch');
const argv = require('minimist')(process.argv.slice(2));

const url = argv.api || 'https://devdata.gramene.org/sorghum_v2';

const idFile = argv.ids;
let targetTaxa = argv.taxa;
let keep = {};
var lineReader = require('readline').createInterface({
  input: require('fs').createReadStream(idFile)
});
// https://data.gramene.org/oryza_v3/search?q={!graph%20from=gene_tree%20to=gene_tree}id:AT2G43490&fq=taxon_id:4558
console.log("id\tortholog.id\tortholog.taxon_id\tortholog.system_name");
lineReader.on('line', function(id) {
  fetch(`${url}/search?rows=10000&fl=id,taxon_id,system_name&q={!graph from=gene_tree to=gene_tree}id:${id}&fq=taxon_id:${targetTaxa}`).then(res => {
    if (res) {
      res.json().then(data => {
        if (data.response) {
          data.response.docs.forEach(doc => {
            console.log([id,doc.id,doc.taxon_id,doc.system_name].join("\t"));
          })          
        }
      })
    }
  })
})
