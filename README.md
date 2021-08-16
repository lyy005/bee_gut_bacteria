# Species formation of gut-restricted bacteria in social bees

This document is a walkthrough of the methods and codes for comparative analysis of the bee gut restricted bactera, *Gilliamella* and *Snodgrassella* genomes. In the paper, we (1) annotated the genomes; (2) performed ortholog assignment and constructed phylogenetic trees; (3) measured gene flow based on PopCOGenT, and (4) performed functional enrichment analysis using anvi'o. Genomes and results of major steps are available on Zenodo: 

## 1 - Gene Annotation

### 1.1 - Annotation using anvi'o db

    # Genomes can be download from Zenodo
    # Install anvio version 7  
    # Make a list of genome names    
    ls *.fna | perl -ne '@a=split/\./; print "$a[0]\n";' > fna.name.list 
    
    # make anvio dbs and gene annotation
    cat fna.name.list | while read i
    do
     echo $i

     echo "step a: rename fasta files"
     anvi-script-reformat-fasta $i\.fna -o $i\-fixed.fa -l 0 --simplify-names
     mv $i\-fixed.fa $i\.fa

     echo "step b: contigs -> database; predict ORF"
     anvi-gen-contigs-database -f $i\.fa -o $i\.db -n "Snodgrassella" # or Gilliamella for Gilliamella genomes

     echo "step c: run HMMs"
     anvi-run-hmms -c $i\.db --num-threads 20

     echo "step d: run NCBI COGs"
     anvi-run-ncbi-cogs -c $i\.db --num-threads 20 --sensitive
     
     echo "step e: export gff files"
     anvi-get-sequences-for-gene-calls -c $i\.db --export-gff3 -o $i\.gff
     
     echo "step f: annotate 16S rRNA genes"
     anvi-get-sequences-for-hmm-hits -c $i\.db --hmm-source Ribosomal_RNAs --gene Bacterial_16S_rRNA -o $i\.16SrRNA
     
     echo "step g: export the protein sequences"
     anvi-get-sequences-for-gene-calls -c $i\.db --get-aa-sequences -o $i\.protein-sequences.fa
     
     echo "step h: annotate gene functions using EGGNOG"
     emapper.py -i $i\.protein-sequences.fa --output $i\.protein-sequences.fa_maNOG --data_dir /PATH_TO_EGGNOG/eggnog-mapper_v2.1.2/ --cpu 32 --override
     
     echo "step i: annotate gene functions using KEGG"
     anvi-run-kegg-kofams -c $i\.db --kegg-data-dir /PATH_TO_KEGG_DB_FROM_ANVIO/KEGG_db -T 12
     
     # KEGG annotations following the tutorial on anvio
     # https://merenlab.org/2018/01/17/importing-ghostkoala-annotations/
     
    done
    
    
### 1.2 - CRISPR annotation using CRISPRCasFinder

    cat cat fna.name.list | while read i | while read i
    do 
    echo $i
    #singularity exec -B $PWD CrisprCasFinder.simg perl /usr/local/CRISPRCasFinder/CRISPRCasFinder.pl -so /usr/local/CRISPRCasFinder/sel392v2.so -cf /usr/local/CRISPRCasFinder/CasFinder-2.0.3 -drpt /usr/local/CRISPRCasFinder/supplementary_files/repeatDirection.tsv -rpts /usr/local/CRISPRCasFinder/supplementary_files/Repeat_List.csv -cas -def G -out ./crispr_out/$i -in $i\.fa -cpuMacSyFinder 8 1> ./crispr_out/$i.1.log 2> ./crispr_out/$i.2.log 
    done

### 1.3 - CAZYme annotation using DBCan2
    # note that the ortholog.all.protein.fasta is from step 2 ortholog assignment
    python run_dbcan.py ortholog.all.protein.fasta protein --out_dir output --dia_cpu 30 --hmm_cpu 30 --tf_cpu 30 --dia_cpu 30 --hotpep_cpu 30

## 2 - Ortholog assignment
    # the create external genome file with path to each genome db (details here: https://merenlab.org/software/anvio/help/main/artifacts/external-genomes/)
    anvi-gen-genomes-storage -e external-genomes.txt -o GGENOMES.db
    
    # ortholog assignment
    anvi-pan-genome -g GGENOMES.db --project-name "G_genomes" --output-dir GGENOMES --num-threads 50 --mcl-inflation 10
    
    # export all annotated coding sequences and protein sequences:
    anvi-get-sequences-for-gene-clusters -p GGENOMES/G_genomes-PAN.db -g GGENOMES.db -o ortholog.all.dna.fasta --report-DNA-sequences --max-num-genes-from-each-genome 1000
    anvi-get-sequences-for-gene-clusters -p GGENOMES/G_genomes-PAN.db -g GGENOMES.db -o ortholog.all.protein.fasta --max-num-genes-from-each-genome 1000
    
    # export ortholog assignment and functions in a table (these files are available on zenodo: XXXX)
    anvi-script-add-default-collection -p GGENOMES/G_genomes-PAN.db
    anvi-summarize -p GGENOMES/G_genomes-PAN.db -g GGENOMES.db -C DEFAULT -o PAN_SUMMARY
    
    

## Citation

Li Y., Leonard S.P., Powell J.E. Moran N.A., 2021. 
