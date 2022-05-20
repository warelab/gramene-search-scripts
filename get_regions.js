#!/usr/bin/env node
var through2 = require('through2');
var byline = require('byline');
var fs = require('fs');
var _ = require('lodash');
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
const cds = argv.cds || 0; // flag for limiting output to cds only
const nofasta = argv.nofasta || 0; // option to turn off fasta output

// read gene list from input
// fetch gene doc from gramene API
// define region to extract
//  If promoter region, use for example: 2kb upstream of first codon of canonical transcript to 500bp downstream
//  If gene region, from start of gene to end of gene + upstream + downstream padding (optional)
// fetch region from ensembl REST API
// if gene region, output gff of canonical transcript with coords adjusted relative to gene start (+ upstream padding)

let gff_writer = fs.createWriteStream(`${outPrefix}.gff`);
gff_writer.write("##gff-version 3\n");
let fasta_writer = fs.createWriteStream(`${outPrefix}.fasta`);

var reader = byline(fs.createReadStream(idFile));

var fetcher = through2.obj(function(line, enc, done) {
  var that = this;
  line = line.toString();
  let cols = line.split("\t");
  const id = cols[0];
  fetch(`${url}/genes?idList=${id}`).then(res => {
    res.json().then(data => {
      if (data[0]) {
        that.push(data[0]);
      }
      done();
    });
  })
});

var doGene = through2.obj(function(gene, enc, done) {
  var that = this;
  let from = gene.location.start;
  let to = gene.location.end;
  let my_up = upstream;
  let my_down = downstream;
  if (mode === "gene") {
    if (upstream > 0) {
      if (gene.location.strand === 1) {
        from -= upstream;
        if (from < 0) {
          my_up += from;
          from = 0;
        }
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
        if (from < 0) {
          my_down += from;
          from = 0;
        }
      }
    }
    // output gff for canonical transcript relative to gene +/- padding
    const gfrom = my_up + 1;
    const gto = gene.location.end - gene.location.start + my_up + 1
    gff_writer.write(`##sequence-region   ${mode.toUpperCase()}.${gene._id} ${gfrom} ${gto}\n`)
    let gff = [
      `${mode.toUpperCase()}.${gene._id}`,
      `panset`,
      `gene`,
      gfrom,
      gto,
      `.`,
      `+`,
      `.`,
      `ID=${gene._id};Name=${gene._id}\n`
    ];
    gff_writer.write(gff.join("\t"));
    const can_id = gene.gene_structure.canonical_transcript;
    const cans = gene.gene_structure.transcripts.filter(t => t.id === can_id);
    gff[2] = `mRNA`;
    if (cds === 1 && cans[0].cds) {
      gff[3] = ggp.remap(gene, cans[0].cds.start, 'transcript', 'gene', can_id) + my_up;
      gff[4] = ggp.remap(gene, cans[0].cds.end, 'transcript', 'gene', can_id) + my_up;
    }
    else {
      gff[3] = ggp.remap(gene, 1, 'transcript', 'gene', can_id) + my_up;
      gff[4] = ggp.remap(gene, cans[0].length, 'transcript', 'gene', can_id) + my_up;
    }
    gff[8] = `ID=${can_id};Name=${can_id};Parent=${gene._id}\n`;
    gff_writer.write(gff.join("\t"));
    // exons
    gff[2] = `exon`;
    let exon_coords = cans[0].exon_junctions || [];
    exon_coords.unshift(0);
    exon_coords.push(cans[0].length);
    exon_coords = _.sortedUniq(exon_coords.sort((a,b) => a-b));
    for(let i=1;i<exon_coords.length;i++) {
      gff[3] = ggp.remap(gene, exon_coords[i-1]+1, 'transcript','gene',can_id) + my_up;
      gff[4] = ggp.remap(gene, exon_coords[i], 'transcript','gene',can_id) + my_up;
      gff[8] = `ID=${can_id}.exon.${i};Parent=${can_id}\n`;
      gff_writer.write(gff.join("\t"));
    }
    if (cans[0].cds) {
      // insert exon junctions for the cds start and end
      exon_coords.push(cans[0].cds.start-1, cans[0].cds.end);
      exon_coords = _.sortedUniq(exon_coords.sort((a,b) => a-b));
      if (cds === 1) {
        exon_coords = exon_coords.filter(ec => ec >= cans[0].cds.start - 1 && ec <= cans[0].cds.end);
      }
      let cds_length_so_far = 0;
      for(let i=1;i<exon_coords.length;i++) {
        gff[2] = 'CDS';
        gff[3] = ggp.remap(gene, exon_coords[i-1]+1, 'transcript','gene',can_id) + my_up;
        gff[4] = ggp.remap(gene, exon_coords[i], 'transcript','gene',can_id) + my_up;
        gff[7] = '.';
        if (exon_coords[i] <= cans[0].cds.start) {
          gff[2] = 'five_prime_UTR';
        }
        else if (exon_coords[i] > cans[0].cds.end) {
          gff[2] = 'three_prime_UTR';
        }
        else {
          gff[7] = cds_length_so_far % 3;
          if (gff[7] > 0) {
            gff[7] = gff[7] === 1 ? 2 : 1;
          }
          cds_length_so_far += gff[4] - gff[3] + 1;
        }
        gff[8] = `ID=${can_id}-${gff[2]};Parent=${can_id}\n`;
        gff_writer.write(gff.join("\t"));
      }
    }
  }
  else {
    
    if (upstream > 0) {
      
    }
  }
  if (nofasta === 0) {
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
          fasta_writer.write(`>${mode.toUpperCase()}.${gene._id} ${gene.system_name} ${gene.location.region}:${from}..${to}:${gene.location.strand}\n${seq}\n`);
          done();
        })
      }
    })  
  }
  else {
    done();
  }
});

reader
.pipe(fetcher)
.pipe(doGene);
