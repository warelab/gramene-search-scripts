#!/usr/bin/env node
var through2 = require('through2');
var byline = require('byline');
var fs = require('fs');
const fetch = require('node-fetch');
const argv = require('minimist')(process.argv.slice(2));
var ggp = require('gramene-gene-positions');
// var genomicStart = ggp.remap(gene, domain.start, 'protein', 'genome');

const idFile = argv.ids;
const outPrefix = argv.prefix;
const url = argv.api || 'https://devdata.gramene.org/sorghum_v2';
const ens = argv.ens || 'https://data.gramene.org/pansite-ensembl';
const mode = argv.mode;
const upstream = argv.upstream || 0;
const downstream = argv.downstream || 0;

// read gene list from input
// fetch gene doc from gramene API
// define region to extract
//  If promoter region, use for example: 2kb upstream of first codon of canonical transcript to 500bp downstream
//  If gene region, from start of gene to end of gene + upstream + downstream padding (optional)
// fetch region from ensembl REST API
// if gene region, output gff of canonical transcript with coords adjusted relative to gene start (+ upstream padding)

let gff_writer = fs.createWriteStream(`${outPrefix}.gff`);
let fasta_writer = fs.createWriteStream(`${outPrefix}.fasta`);

var reader = byline(fs.createReadStream(idFile));

var fetcher = through2.obj(function(line, enc, done) {
  var that = this;
  line = line.toString();
  let cols = line.split("\t");
  const id = cols[0];
  fetch(`${url}/genes?idList=${id}`).then(res => {
    res.json().then(data => {
      that.push(data[0]);
      done();
    });
  })
});

var doGene = through2.obj(function(gene, enc, done) {
  var that = this;
  let from = gene.location.start;
  let to = gene.location.end;
  if (mode === "gene") {
    if (upstream > 0) {
      if (gene.location.strand === 1) {
        from -= upstream;
      }
      else {
        to += upstream;
      }
    }
    if (downstream > 0) {
      if (gene.location.strand === 1) {
        to += downstream;
      }
      else {
        from -= downstream;
      }
    }
    // output gff for canonical transcript relative to gene +/- padding
    let gff = [
      `${mode.toUpperCase()}.${gene._id}`,
      `panset`,
      `gene`,
      upstream + 1,
      gene.location.end - gene.location.start + upstream + 1,
      `.`,
      `+`,
      `.`,
      `ID=${gene._id};Name=${gene._id}\n`
    ];
    gff_writer.write(gff.join("\t"));
    const can_id = gene.gene_structure.canonical_transcript;
    const cans = gene.gene_structure.transcripts.filter(t => t.id === can_id);
    gff[2] = `mRNA`;
    gff[3] = ggp.remap(gene, 1, 'transcript', 'gene', can_id) + upstream;
    gff[4] = ggp.remap(gene, cans[0].length, 'transcript', 'gene', can_id) + upstream;
    gff[8] = `ID=${can_id};Name=${can_id};Parent=${gene._id}\n`;
    gff_writer.write(gff.join("\t"));
    // exons
    gff[2] = `exon`;
    let exon_coords = cans[0].exon_junctions || [];
    exon_coords.unshift(0);
    exon_coords.push(cans[0].length);
    for(let i=1;i<exon_coords.length;i++) {
      gff[3] = ggp.remap(gene, exon_coords[i-1]+1, 'transcript','gene',can_id) + upstream;
      gff[4] = ggp.remap(gene, exon_coords[i], 'transcript','gene',can_id) + upstream;
      gff[8] = `ID=${can_id}.${i};Parent=${can_id}\n`;
      gff_writer.write(gff.join("\t"));
    }
  }
  else {
    
    if (upstream > 0) {
      
    }
  }
  let re = /\+/g;
  gene.location.region = gene.location.region.replace(re, "p");
  const ensURL = `${ens}/sequence/region/${gene.system_name}/${gene.location.region}:${from}..${to}:${gene.location.strand}?content-type=application/json`;
  fetch(ensURL).then(ens_res => {
    if (ens_res) {
      ens_res.json().then(seq_data => {
        let seq = seq_data.seq;
        if (!seq) {
          throw(`seq not defined for for ${gene._id} in response from ${ensURL}`)
        }
        if (gene.location.strand === -1) {
          seq = seq_data.seq.split("").reverse().join("");
        }
        fasta_writer.write(`>${mode.toUpperCase()}.${gene._id} ${gene.system_name} ${gene.location.region}:${from}..${to}:${gene.location.strand}\n${seq}\n`);
        done();
      })
    }
  })  
});

reader
.pipe(fetcher)
.pipe(doGene);

// lineReader.on('line', function(id) {
//   fetch(`${url}/genes?idList=${id}`).then(res => {
//     if (res) {
//       res.json().then(data => {
//         let gene = data[0];
//         let from = gene.location.start;
//         let to = gene.location.end;
//         if (mode === "gene") {
//           if (upstream > 0) {
//             if (gene.location.strand === 1) {
//               from -= upstream;
//             }
//             else {
//               to += upstream;
//             }
//           }
//           if (downstream > 0) {
//             if (gene.location.strand === 1) {
//               to += downstream;
//             }
//             else {
//               from -= downstream;
//             }
//           }
//           // output gff for canonical transcript relative to gene +/- padding
//           let gff = [
//             `${mode.toUpperCase()}:${gene._id}`,
//             `panset`,
//             `gene`,
//             upstream + 1,
//             gene.location.end - gene.location.start + upstream + 1,
//             `.`,
//             `+`,
//             `.`,
//             `ID=${gene._id};Name=${gene._id}`
//           ];
//           console.log(gff.join("\t"));
//           const can_id = gene.gene_structure.canonical_transcript;
//           const cans = gene.gene_structure.transcripts.filter(t => t.id === can_id);
//           gff[2] = `mRNA`;
//           gff[3] = ggp.remap(gene, 1, 'transcript', 'gene', can_id) + upstream;
//           gff[4] = ggp.remap(gene, cans[0].length, 'transcript', 'gene', can_id) + upstream;
//           gff[8] = `ID=${can_id};Name=${can_id};Parent=${gene._id}`;
//           console.log(gff.join("\t"));
//           // exons
//           gff[2] = `exon`;
//           let exon_coords = cans[0].exon_junctions || [];
//           exon_coords.unshift(0);
//           exon_coords.push(cans[0].length);
//           for(let i=1;i<exon_coords.length;i++) {
//             gff[3] = ggp.remap(gene, exon_coords[i-1]+1, 'transcript','gene',can_id) + upstream;
//             gff[4] = ggp.remap(gene, exon_coords[i], 'transcript','gene',can_id) + upstream;
//             gff[8] = `ID=${can_id}.${i};Parent=${can_id}`;
//             console.log(gff.join("\t"));
//           }
//         }
//         else {
//
//           if (upstream > 0) {
//
//           }
//         }
//         const ensURL = `${ens}/sequence/region/${gene.system_name}/${gene.location.region}:${from}..${to}:${gene.location.strand}?content-type=application/json`;
//         // console.log(ensURL);
//         fetch(ensURL).then(ens_res => {
//           if (ens_res) {
//             ens_res.json().then(seq_data => {
//               let seq = seq_data.seq;
//               if (gene.location.strand === -1) {
//                 seq = seq_data.seq.split("").reverse().join("");
//               }
//               console.log(`>${mode.toUpperCase()}:${gene._id}\n${seq}`);
//             })
//           }
//         })
//       })
//     }
//   })
// })
//
