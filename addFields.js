#!/usr/bin/env node

var fs = require('fs');
var _ = require('lodash');
var searchInterface = require('./index.js').client;

// addFields.js input.txt idColumn collection {'field':'_id', 'fl':['name','species'], 'taxon_id':3702}
/*
{
  field: the field containing the query terms
  fl: list of fields to return (excluding query field)
  <any other key-value pairs>
}
*/
var inputFile = process.argv[2];
if (inputFile === '-') inputFile = '/dev/fd/0';
var idColumn = process.argv[3];
var collection = process.argv[4];
var etc = process.argv.slice(5).join(' ');
var params = etc ? JSON.parse(etc) : {};
var fields = params.fl;
var inputData = [];
var idSet={};
require('readline').createInterface(
  {
    input: fs.createReadStream(inputFile),
    terminal: false
  }
).on('line', function(line)
  {
    var cols = line.split('\t');
    var id = cols[idColumn];
    inputData.push(cols);
    idSet[id]=1;
  }
).on('close', function()
  {
    var idList = Object.keys(idSet);
    searchInterface.coreLookup(collection,idList,params).then(function(lut) {
      inputData.forEach(function(cols) {
        var id = cols[idColumn];
        fields.forEach(function(f) {
          var v = '-';
          if (lut[id]) {
            v = lut[id].map(function(doc) {
              return doc[f];
            });
          }
          cols.push(_.uniq(v));
        });
        console.log(cols.join("\t"));
      });
    });
  }
);
