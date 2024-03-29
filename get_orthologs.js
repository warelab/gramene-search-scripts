#!/usr/bin/env node
const fetch = require('node-fetch');
const argv = require('minimist')(process.argv.slice(2));

const url = argv.api || 'https://devdata.gramene.org/sorghum_v2';

const idFile = argv.ids;
const targetTaxa = argv.taxa;
let keep = {};
if (targetTaxa) {
  targetTaxa.split(',').forEach(tid => keep[tid]=1)
}
var lineReader = require('readline').createInterface({
  input: require('fs').createReadStream(idFile)
});

console.log("id\tortholog.id\tortholog.taxon_id\tortholog.system_name");
lineReader.on('line', function(id) {
  fetch(`${url}/search?rows=10000&fl=id,taxon_id,system_name&q=homology__all_orthologs:${id}`).then(res => {
    res.json().then(data => {
      data.response.docs.forEach(doc => {
        if (!targetTaxa || keep[doc.taxon_id]) {
          console.log([id,doc.id,doc.taxon_id,doc.system_name].join("\t"));
        }
      })
    })
  })
})
