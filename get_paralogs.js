#!/usr/bin/env node
const fetch = require('node-fetch');
const argv = require('minimist')(process.argv.slice(2));

const url = argv.api || 'https://devdata.gramene.org/sorghum_v2';

const idFile = argv.ids;

var lineReader = require('readline').createInterface({
  input: require('fs').createReadStream(idFile)
});

console.log("id\tparalog.id\ttaxon_id");
lineReader.on('line', function(id) {
  fetch(`${url}/search?rows=10000&fl=id,taxon_id&q=homology__within_species_paralog:${id}`).then(res => {
    res.json().then(data => {
      data.response.docs.forEach(doc => {
        console.log([id,doc.id,doc.taxon_id].join("\t"));
      })
    })
  })
})
