#parse genbank file into format that can be read by database (genes table)

from Bio import SeqIO
import gzip

geneCount = 0
exonCount = 0
repCount = 0
replicontbl = []
genetbl = []
exontbl = []


files = ["GCF_000005845.2_ASM584v2_genomic.gbff.gz","GCF_000576515.1_ASM57651v1_genomic.gbff.gz"]


#Open each genome in the file list
for genome_ID, f in enumerate(files):
  with gzip.open(f,'rb') as inzip:

    #parse each record as a replicon, store in replicon list
    for record in SeqIO.parse(inzip, "genbank"):
      repentry = []

      #replicon ID
      repCount += 1
      repentry.append(repCount)

      #genome_ID
      repentry.append(genome_ID+1)

      #name
      repentry.append(record.description)

      #accession
      repentry.append(record.id)

      #rep_type
      if "chromosome" in record.description.lower():
        repentry.append("chromosome")
      else:
        repentry.append("plasmid")

      #rep_struct
      repentry.append(record.annotations['topology'])

      replicontbl.append(repentry)


      #parse each feature that is a CDS, store each as a list in gene table
      for feat in record.features:
        if feat.type == "CDS":

          try:
            int(feat.location.start)
            int(feat.location.end)
          except ValueError:
            continue;

          entry = []

          #gene_id
          geneCount += 1
          entry.append(geneCount)
          
          #genome_id
          entry.append(genome_ID+1)

          #replicon_id
          entry.append(repCount)

          #locus_tag
          entry.append(feat.qualifiers['locus_tag'][0])

          #gene_name
          if "gene" in feat.qualifiers:
            entry.append(feat.qualifiers['gene'][0])
          else:
            entry.append("NULL")

          #strand
          if feat.location.strand > 0:
            entry.append("+")
          else:
            entry.append("-")

          #num_exons
          entry.append(len(feat.location.parts))

          #length
          entry.append(len(feat.location))

          #product_name
          if "product" in feat.qualifiers:
            entry.append(feat.qualifiers['product'][0])
          else:
            entry.append("NULL")

          genetbl.append(entry)

          #parse location(s) for each feature as an exon, insert into exon list
          for loc in feat.location.parts:
            
            exentry = []
            exonCount += 1

            #gene_id
            exentry.append(geneCount)
            #exon
            exentry.append(exonCount)
            #left_pos
            exentry.append(loc.start)
            #right_pos
            exentry.append(loc.end)
            #length
            exentry.append(len(loc))

            exontbl.append(exentry)



#write each list to a separate tab delimited file for loading into database
with open("genesout.txt",'w') as out:
  for entry in genetbl:
    out.write('\t'.join([str(x) for x in entry])+'\n')

with open("exonsout.txt",'w') as out:
  for entry in exontbl:
    out.write('\t'.join([str(x) for x in entry])+'\n')

with open("repsout.txt",'w') as out:
  for entry in replicontbl:
    out.write('\t'.join([str(x) for x in entry])+'\n')

print "success"
print str(exonCount) + " exons"
print str(geneCount) + " genes"
print str(repCount) + " replicons"