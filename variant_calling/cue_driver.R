f = "/stanley/levin_asap_storage/kwanho/Revio_gDNA/bam/filelist.tsv"
tab = read.table(f, sep='\t', header=F)
colnames(tab) = c("sample", "bam")

# check total size before downloading
#all_bams = paste0("gcloud storage ls -L ", paste(tab$bam, collapse=' '), " | grep "[0-9]GiB" | awk '{print $6}' | awk '{gsub(/[()GiB]/, "", $1); sum+=$1} END {print sum}' > total_bam_size.txt")
#write.table(all_bams, "check_size_all_bams.sh", sep='\t', quote=F, row.names=F, col.names=F)

myconf = data.frame(matrix(, nrow=5, ncol=1, 
	dimnames=list(c('bam', 'chr_names', 'model_path','fa','n_cpus'),c('value'))))
myconf['chr_names',] = "null"
myconf['model_path',] = '"/stanley/levin_asap_storage/kwanho/Revio_gDNA/variant_calls/cue2/models/cue2.hifi.v1.pt"'
myconf['fa',] = '"/stanley/levin_asap_storage/kwanho/Revio_gDNA/refs/no_alt_GRCh38_ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"'
myconf['n_cpus',] = 16
myconf.base = myconf

run_script = "/stanley/levin_asap_storage/kwanho/Revio_gDNA/variant_calls/cue2/run_cue_template.sh"

# configure for Cue2 run
for (i in 32:nrow(tab)) {
	samp = tab[i, 'sample']
	bam_gcs = tab[i, 'bam']
	bai_gcs = paste0(bam_gcs, ".bai")
	print(samp)
	
	dir.create(samp, showWarnings=F)
	setwd(samp)

	out_bam = paste0(samp, ".bam")
	out_bai = paste0(out_bam, ".bai")
	
	# Get BAM from GCP
	system(paste("gcloud storage cp", bam_gcs, out_bam))
	system(paste("gcloud storage cp", bai_gcs, out_bai))

	# config file
	myconf = myconf.base
	myconf['bam',] = paste0('"', out_bam ,'"')
	write.table(myconf, "config.yaml", sep=': ', quote=F, row.names=T, col.names=F)

	# copy driver script
	system(paste("cp", run_script, "driver.sh"))

	# submit job
	system(paste0("qsub -N cue_", samp, " -j y driver.sh"))

	setwd('../')
}

