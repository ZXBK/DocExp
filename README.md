# DocExp
Execution of DEgene calculation

Input: Genome fasta file & Genome gff file & RNA raw seqeuences downloaded from NCBI

Output: FPKM table

## Docker hub: https://hub.docker.com/r/lsbnb/docexpress

### Usage:
docker pull lsbnb/docexpress

docker run -itd --name docexp -p 8999:80 -p 8021:21 -p 8022:22 -v Path:/in -w /in lsbnb/docexpress bash

docker exec -it docexp bash

chmod -R 777 /in/

nohup sh doc.sh -p 20 -d /path/to/sra/ -g genome.fna -f genome.gff &
