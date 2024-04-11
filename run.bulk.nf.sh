
#(base) biolab@BioLab:/mnt/second/yasinkaymaz$ 
nextflow run nf-core/rnaseq -profile docker \
--input samplesheet.2.csv \
--outdir results.1 \
--max_cpus 90 \
--max_memory 300.GB \
--fasta /home/sharedFolder/humanSTARindex/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--gtf /home/sharedFolder/humanSTARindex/Homo_sapiens.GRCh38.99.gtf \
--star_index /home/sharedFolder/humanSTARindex/ \
-resume \
--skip_qualimap \
--skip_rseqc \
--skip_deseq2_qc