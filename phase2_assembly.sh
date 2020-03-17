#!/bin/bash
#1: Sample Name
#2: Working Directory
#3: Threads
#4: Illumina 1
#5: Illumina 2
#6: Nanopore
#7: Genome Size


  #### Pilon_func
  # Arguments:
  #  1. input
  #  2. working_dir
  #  3. log
  #  4. pilon_mode
  #  5. illumina_r1_fastq
  #  6. illumina_r2_fastq
  #  7. threads
  pilon_func() {
    if [ $# != 7 ]
    then
      echo 'Command should be run as:'
      echo 'pilon_func {Assembly} {Working Directory} {Log File} {Pilon Mode: (paired|frag)} {Illumina R1} {Illumina R2} {Threads}'
      return 1
    fi
    local current_asm_fa=$1
    local working_dir=$2
    local log=$3
    local pilon_mode=$4
    local illumina_r1_fastq=$5
    local illumina_r2_fastq=$6
    local threads=$7

    echo  "$(date +%T%Z) - Pilon: Starting" >> $log

    # pilon polishing
    if [ "$pilon_mode" != 'None' ]; then
        local illumina_map="$working_dir/illumina_to_asm.bam"
        local new_asm_fa="$working_dir/pilon_finished.assembly.fasta"  # .fasta is required pilon suffix
        bwa index "$current_asm_fa"
        if [ "paired" = "$pilon_mode" ]; then
            bwa mem -t $threads "$current_asm_fa" "$illumina_r1_fastq" "$illumina_r2_fastq" | samtools view -b - > "$illumina_map"
            samtools sort -o "$illumina_map".sorted "$illumina_map" && mv "$illumina_map".sorted "$illumina_map"
            samtools index "$illumina_map"
            pilon -Xmx120g --genome "$current_asm_fa" --jumps "$illumina_map" --changes --verbose --nostrays --fix all --output "${new_asm_fa%.fasta}"
        elif [ "frag" = "$pilon_mode" ]; then
            bwa mem "$current_asm_fa" "$illumina_r1_fastq" | samtools view -b - > "$illumina_map"
            pilon --genome "$current_asm_fa" --frags "$illumina_map" --output "${new_asm_fa%.fasta}"
        fi
        current_asm_fa=$new_asm_fa
    fi

    echo  "$(date +%T%Z) - Pilon: Finished" >> $log

    pilon_output="$current_asm_fa"
    return 0
  }

########################## Main Pipeline #####################################
if [ $# != 7 ]
then
  echo 'Command should be run as:'
  echo 'phase2_assembly.sh {Sample Name} {Working Directory} {Threads} {Illumina R1 Path} {Illumina R2 Path} {ONT Path} {Expected Genome Size}'
  echo 'The script will create a directory named {Sample Name} in the {Working Directory}'
else
  sample="${1}"
  wd="${2}/${sample}"
  threads=${3}
  ill1="${4}"
  ill2="${5}"
  ONT="${6}"
  gsize="${7}"

  # Spades assembly
  ill_asm="${wd}/ASM_work/spades/scaffolds.fasta"
  # Polished flye assembly
  ONT_asm="${wd}/ASM_work/flye/pilon09/pilon_finished.assembly.fasta"
  # Non-polished flye assembly
  flye_asm="${wd}/ASM_work/flye/assembly.fasta"
  # Blast results of spades vs flye assemblies
  blast_xml="${wd}/ASM_work/blast_results.xml"
  # Contigs in spades assembly not contained in flye assembly
  spades_only="${wd}/ASM_work/spades_only.fasta"
  # Map of Oxford reads against spades_only
  minimap_paf="${wd}/ASM_work/spades_only_ont.paf"
  # Ids of Spades_only contigs that are contained in Oxford reads
  spades_mapped_ids="${wd}/ASM_work/spades_mapped_ids.txt"
  # Ids of Oxford reads containing spades contigs
  ONT_spades_ids="${wd}/ASM_work/ONT_spades_ids.txt"
  # Spades_only contigs that are not contained in Oxford reads
  spades_unmapped_contigs="${wd}/ASM_work/spades_unmapped_contigs.fasta"
  # Oxford reads that contain Spades_only contigs
  ONT_context_contigs="${wd}/ASM_work/ONT_context_contigs.fasta"
  # Polished ONT_context_contigs
  polished_context_contigs="${wd}/ASM_work/pilon09/pilon_finished.assembly.fasta"
  # Final assembly containing flye assembly, spades contigs not contained by flye or Oxford reads and polished Oxford reads containing spades contigs
  final_asm="${wd}/ASM_work/${sample}.fasta"

  # Setup directory structure
  mkdir $wd
  # mkdir ${wd}/ONT/
  # mkdir ${wd}/ONT/gDNA
  # mkdir ${wd}/ONT/RNA
  # mkdir ${wd}/Ill/
  # mkdir ${wd}/Ill/gDNA
  # mkdir ${wd}/Ill/RNA

  # Copy input files into working directory
  # if [ ! -f ${wd}/Ill/gDNA/${sample}_R1.fastq.gz ] || [ ! -f ${wd}/Ill/gDNA/${sample}_R2.fastq.gz ];
  # then
  #   cp ${ill1} ${wd}/Ill/gDNA/${sample}_R1.fastq.gz
  #   cp ${ill2} ${wd}/Ill/gDNA/${sample}_R2.fastq.gz
  # fi
  # if [ ! -f ${wd}/ONT/gDNA/${sample}_ONT.fastq.gz ];
  # then
  #   cp ${ONT} ${wd}/ONT/gDNA/${sample}_ONT.fastq.gz
  # fi

  # Set input variables to the new paths
  # ill1=${wd}/Ill/gDNA/${sample}_R1.fastq.gz
  # ill2=${wd}/Ill/gDNA/${sample}_R2.fastq.gz
  # ONT=${wd}/ONT/gDNA/${sample}_ONT.fastq.gz

  # Setup directories for assembly and combining things
  mkdir ${wd}/ASM_work
  mkdir ${wd}/ASM_work/spades
  mkdir ${wd}/ASM_work/flye

  # Run spades assembly
  if [ ! -f ${ill_asm} ];
  then
    spades.py -1 ${ill1} -2 ${ill2} --nanopore ${ONT} -t ${threads} -o ${wd}/ASM_work/spades --meta
  fi

  # Run flye assembly
  if [ ! -f ${flye_asm} ];
  then
    flye --nano-raw ${ONT} -g ${gsize} -t ${threads} --plasmids --meta  -o ${wd}/ASM_work/flye
  fi

  # Polish flye assembly 10 times
  if [ ! -f ${ONT_asm} ];
  then
    asm=${wd}/ASM_work/flye/assembly.fasta
    for r in {00..09};
    do
      mkdir ${wd}/ASM_work/flye/pilon${r}
      pilon_func ${asm} ${wd}/ASM_work/flye/pilon${r} ${wd}/ASM_work/flye/pilon.log paired ${ill1} ${ill2} ${threads}
      asm=${wd}/ASM_work/flye/pilon${r}/pilon_finished.assembly.fasta
    done
  fi

  # Blast spades contigs against flye assembly
  if [ ! -f ${blast_xml} ];
  then
    blastn -query ${ill_asm} -subject ${ONT_asm} -out ${blast_xml} -outfmt "6 std qlen slen"
  fi

  # Remove spades contigs contained by flye contigs
  if [ ! -f ${spades_only} ];
  then
    get_uncovered_contigs.py -b ${blast_xml} -q ${ill_asm} -o ${spades_only}
  fi

  # Map spades only contigs to Oxford reads
  if [ ! -f ${minimap_paf} ];
  then
    minimap2 --paf-no-hit -t ${threads} -x map-ont -N 1 ${ONT} ${spades_only} > ${minimap_paf}
  fi

  # Get spades only contig ids that are covered by Oxford reads
  if [ ! -f ${spades_mapped_ids} ];
  then
    cat ${minimap_paf} | cut -f1,2,6,11 | awk '{if ($4/$2 >= 0.9){print $0}}' |  cut -f1 | sort | uniq > ${spades_mapped_ids}
  fi

  # Get Oxford reads' ids that contain spades' ids
  if [ ! -f ${ONT_spades_ids} ];
  then
     cat ${minimap_paf} | cut -f1,2,6,11 | awk '{if ($4/$2 >= 0.9){print $0}}' |  cut -f3 | sort | uniq > ${ONT_spades_ids}
  fi

  # Get Oxford reads that contain spades contigs and spades contigs that are not contained by Oxford reads
  if [ ! -f ${spades_unmapped_contigs} ] || [ ! -f ${ONT_context_contigs} ];
  then
    get_spades_context.py --spades_ids ${spades_mapped_ids} --ONT_ids ${ONT_spades_ids} \
    --spades_file ${spades_only} --ONT_file ${ONT} --spades_output ${spades_unmapped_contigs} --ont_output ${ONT_context_contigs}
  fi

  # Polish Oxford reads containing spades contigs
  if [ ! -f ${polished_context_contigs} ];
  then
    asm=${ONT_context_contigs}
    for r in {00..09};
    do
      mkdir ${wd}/ASM_work/pilon${r}
      pilon_func ${asm} ${wd}/ASM_work/pilon${r} ${wd}/ASM_work/pilon.log paired ${ill1} ${ill2} ${threads}
      asm=${wd}/ASM_work/pilon${r}/pilon_finished.assembly.fasta
    done
  fi

  # Combine polished Oxford reads, flye contigs, and spades contigs not covered by flye nor by Oxford reads
  if [ ! -f ${final_asm} ];
  then
    combine_fasta.py --sample ${sample} --unmapped ${spades_unmapped_contigs} --context ${polished_context_contigs} --flye ${ONT_asm} --output ${final_asm}
  fi


fi
