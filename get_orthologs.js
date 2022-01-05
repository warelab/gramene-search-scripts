#!/usr/bin/env node
const fetch = require('node-fetch');
const nReadlines = require('n-readlines');
const argv = require('minimist')(process.argv.slice(2));

const url = argv.api || 'https://devdata.gramene.org/sorghum_v2';

const idFile = argv.ids;

const idLines = new nReadlines(idFile);

console.log("id\tortholog.id\tortholog.taxon_id\tortholog.system_name");
while (line = idLines.next()) {
  const id = line.toString('ascii');
  fetch(`${url}/search?rows=10000&fl=id,taxon_id,system_name&q=homology__all_orthologs:${id}`).then(res => res.json().then(data => {
    data.response.docs.forEach(doc => {
      console.log([id,doc.id,doc.taxon_id,doc.system_name].join("\t"));
    })
  }));    
}
