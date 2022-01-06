#!/usr/bin/env node
const fetch = require('node-fetch');
const argv = require('minimist')(process.argv.slice(2));

const url = argv.url || 'https://data.gramene.org/vitis1/search';
const facet = argv.facet || 'gene_tree';
const params = {
  'q': argv.q || 'domains__ancestors:2182',
  'facet.field': `{!facet.limit=10000 facet.mincount=1 key=${facet}}${facet}`,
  'rows':0
};
function queryString(params) {
  let qps = [];
  for (const field in params) {
    qps.push(`${field}=${params[field]}`)
  }
  return qps.join("&");
}
fetch(`${url}?${queryString(params)}`, {
  "headers": {
    "accept": "*/*",
    "accept-language": "en-US,en;q=0.9",
    "sec-ch-ua": "\" Not A;Brand\";v=\"99\", \"Chromium\";v=\"96\", \"Google Chrome\";v=\"96\"",
    "sec-ch-ua-mobile": "?0",
    "sec-ch-ua-platform": "\"macOS\"",
    "sec-fetch-dest": "empty",
    "sec-fetch-mode": "cors",
    "sec-fetch-site": "same-site",
    "Referer": "https://vitis.gramene.org/",
    "Referrer-Policy": "strict-origin-when-cross-origin"
  },
  "body": null,
  "method": "GET"
})
.then(res => res.json().then(data => {
  const trees = data.facet_counts.facet_fields[facet].filter((element, index) => index % 2 === 0);
  console.log(trees.join("\n"));
}));
