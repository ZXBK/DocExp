#!/bin/bash
#Kao Chih Hsin
#Feb,2022
###docker pull lsbnb/docexpress_fastqc
###docker run -itd -p 8999:80 -p 8021:21 -p 8022:22 -v path:/in -w /in
### must map the path to /in !!!
###chmod -R 777 /in/
PATH=/root/galaxy/database/dependencies/_conda/envs/__trim-galore@0.4.3/bin/:/root/miniconda2/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin ; export PATH

usage () {
    cat <<HELP_USAGE
    Usage: $0  [-p threads] [-d path ] [-g genome.fna] [-f genome.gff] 
    -p:  threads
    -d:  path ex: /in/sra (path must contained RNA paired file ex. SRRxxxx_1.fastq  )
    -g:  genome file (.fna)
    -f: gff file 
HELP_USAGE
}

### get options
while getopts 'p:d:g:f:h:' OPT; do
        case $OPT in
                p)
                        threads="$OPTARG";;
                d)
                        path="$OPTARG";;
                g)
                        genome="$OPTARG";;
                f)
                        gff="$OPTARG";;
                h)
                        usage;;
                *)
                        usage
                        exit 1;;
        esac
done

### Test if options are provided
if [ -z $threads ];then
    echo  "Option -p is empty"
    usage
    exit 1
fi
if [ -z $path ];then
    echo  "Option -d is empty"
    usage
    exit 1
fi
if [ -z $genome ];then
    echo  "Option -g is empty"
    usage
    exit 1
fi
if [ -z $gff ];then
    echo  "Option -f is empty"
    usage
    exit 1
fi

### Test if files are exist/are provided
#-s FILE exists and has a size greater than zero.
if [ ! -s $genome ];then
    echo -e "The candidate file: '$genome' didn't exists or is empty\nExit\n"
    usage
    #exit 1
fi
if [ ! -s $gff ];then
    echo -e "The candidate file: '$gff' didn't exists or is empty\nExit\n"
    usage
    exit 1
fi
if [ ! -d ${path} ] 
then
    echo "path didn't exists or is empty"
    usage
    exit 1 
fi

### check TrimGalore program
if [ ! -d /in/TrimGalore-0.6.6/ ];then 
    curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.tar.gz -o /in/trim_galore.tar.gz
    tar xvzf /in/trim_galore.tar.gz
fi

echo using threads: $threads
echo sra path: $path
echo genome: $genome
echo gff file: $gff 
echo 

### executing trimgalore ###
echo ========== Start TrimGalore v0.4.3.1 ==========
echo $(date '+%d %B %T')

if [ ! -d /in/trimG ] ; then mkdir /in/trimG ; fi
cd /in/
for ID in $(ls ${path} | grep -P -o '.*(?=_1.fastq)')
do
    if [ ! -f /in/trimG/${ID}_1.fastq_trimming_report.txt ] ; then 
        nohup /in/TrimGalore-0.6.6/trim_galore -j ${threads} -o /in/trimG/  --paired ${path}${ID}_1.fastq ${path}${ID}_2.fastq  &
        pids="$!"
        echo
        echo $(date '+%d %B %T') ...... Processing TrimGalore PID: $pids
        wait $pids
    fi 
done

sra_file=$(ls ${path} | grep -P -o '.*(?=.fastq)' | wc | awk '{print $1}')
echo "##### sra file number :  ${sra_file} (r1+r2) #####"
trimG_file=$(ls /in/trimG/ | wc | awk '{print $1}')
echo "##### trimG output file number :  ${trimG_file} (r1+r2+report_r1_report_r2) #####"
if [ ${trimG_file} != $((${sra_file}*2)) ] ; then 
    echo " trimG file number wrong "
    exit 1
fi

### executing hisat2 v2.1 ###
echo ========== Start Hisat2 v2.1 ==========
cd /in/
if [ ! -d /in/hisat_build/ ] ; then mkdir /in/hisat_build ; fi
if [ ! -f /in/hisat_build/Species_Index.1.ht* ] ;then
    echo starting hisat2 build process 
    nohup hisat2-build ${genome} /in/hisat_build/Species_Index  &
    pids="$!"
    echo
    echo $(date '+%d %B %T') ...... Processing Hisat2 building PID: $pids
    wait $pids
fi 
if [ ! -f /in/hisat_build/Species_Index.1.ht* ] ; then 
    echo " hisat2 index wrong "
    exit 1 
fi 

if [ ! -d /in/hisat_result ] ; then mkdir /in/hisat_result ; fi
for ID in $(ls ${path} | grep -P -o '.*(?=_1.fastq)')
do
    if [ ! -f /in/hisat_result/${ID}.sam ] ; then 
        nohup hisat2 -p ${threads} -x /in/hisat_build/Species_Index -1 /in/trimG/${ID}_1_val_1.fq -2 /in/trimG/${ID}_2_val_2.fq \
        -S /in/hisat_result/${ID}.sam  &
        pids="$!"
        echo
        echo $(date '+%d %B %T') ...... Processing Hisat2 PID: $pids
        wait $pids
        rm /in/trimG/${ID}_1_val_1.fq /in/trimG/${ID}_2_val_2.fq
    fi
done

sam_file=$(ls /in/hisat_result/ | wc | awk '{print $1}')
echo "##### sam file number :  ${sam_file}  #####"
if [ ${sam_file} != `expr ${sra_file} / 2 `  ] ; then
    echo " sam file number wrong "
    exit 1
fi

for ID in $(ls ${path} | grep -P -o '.*(?=_1.fastq)')
do
    if [ ! -f /in/hisat_result/${ID}_sort.bam ] ; then 
        nohup samtools view -@ ${threads} -b -S /in/hisat_result/${ID}.sam > /in/hisat_result/${ID}.bam 2>> log_samtoolsView &
        pids="$!"
        echo
        echo $(date '+%d %B %T') ...... Processing samtools view ${ID} PID: $pids
        wait $pids
        rm /in/hisat_result/${ID}.sam
        nohup samtools sort -@ ${threads} /in/hisat_result/${ID}.bam /in/hisat_result/${ID}_sort & 
        pids="$!"
        echo
        echo $(date '+%d %B %T') ...... Processing samtools sort ${ID} PID: $pids
        wait $pids
        rm /in/hisat_result/${ID}.bam
    fi 
done


### executing stringtie v1.3.3 ###
echo ========== Start stringtie v1.3.3 ==========

if [ ! -d /in/stringtie_result ] ; then mkdir /in/stringtie_result ; fi
for ID in $(ls ${path} | grep -P -o '.*(?=_1.fastq)')
do
    if [ ! -f /in/stringtie_result/${ID}.table ] ; then 
        nohup stringtie /in/hisat_result/${ID}_sort.bam -G ${gff} \
        -o /in/stringtie_result/${ID}.gtf -p ${threads} -A /in/stringtie_result/${ID}.fpkm -e  &
        pids="$!"
        echo
        echo $(date '+%d %B %T') ...... Processing stringtie PID: $pids
        wait $pids
        rm -rf /in/stringtie_result/tmp*
    fi
done

### merge fpkm table ###
echo ========== Start merging by using csvtk ==========
if [ ! -f /in/csvtk_linux_amd64.tar.gz ] ; then 
    curl -fsSL https://github.com/shenwei356/csvtk/releases/download/v0.24.0/csvtk_linux_amd64.tar.gz -o /in/csvtk_linux_amd64.tar.gz && \
    tar zxvf /in/csvtk_linux_amd64.tar.gz && cp /in/csvtk /usr/local/bin/
fi

if [ ! -d /in/merge_table ] ; then mkdir /in/merge_table ; fi
for ID in $(ls /in/stringtie_result/ | grep -P -o '.*(?=.fpkm)')
do
    if [ ! -f /in/merge_table/${ID}.table ] ; then 
        sort -k1,1 /in/stringtie_result/${ID}.fpkm > /in/merge_table/${ID}.sortfpkm && \
        awk '{print $1,$8}' /in/merge_table/${ID}.sortfpkm > /in/merge_table/${ID}.table
        sed -i 's/ /,/g' /in/merge_table/${ID}.table
        sed -i 's/Gene/KEY/g' /in/merge_table/${ID}.table
        sed -i "s/End/${ID}/g" /in/merge_table/${ID}.table
    fi
done

/in/csvtk join -f 1 /in/merge_table/*.table -O --na 0 > /in/merge_table/fpkm.table

echo $(date '+%d %B %T')
echo ========== Finish :D ==========
