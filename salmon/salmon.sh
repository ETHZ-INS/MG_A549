#! /bin/bash

ref=/reference/Homo_sapiens/Ensembl/GRCh38.p13/Annotation/Release_105-2021
for f in ../*/raw_data/*/; do
  f2=`basename $f | sed 's/pl_//'`;
  echo $f2;
  /common/salmon-1.7.0/bin/salmon quant -p 10 -i $ref/salmon1.7.0/ \
    -l A -1 "$f"*_1.fq.gz -2 "$f"*_2.fq.gz --validateMappings \
    -q -g $ref/tx2gene -o $f2
done

