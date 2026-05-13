nextflow.enable.dsl=2
params.timestamp = ""
params.qc_type = "dna" // Options: "dna" or "rna"

process QC_PLOTS {   
    tag "qc_plots"
	publishDir "${publish_dir}_${params.timestamp}/tertiary_analyses/qc_plots", enabled:"$enable_publish"

	input:
		path rds_file
		path seg_copy_file
		path metrics_file
		path selected_metrics
		path metadata_file
		path plot_qc_config
		val(publish_dir)
		val(enable_publish)

	output:
		path "AllSample-GinkgoSegmentSummary.txt"
		path "nf-preseq-pipeline_all_metrics_mqc.txt", emit : allmetrics_with_cnv
		path "nf-preseq-pipeline_selected_metrics_mqc.txt", emit : selectedmetrics_with_cnv, optional: true
		path "DNA-QC_ConsensusScores.txt", emit : consensus_scores
		path "DNA-QC_ConsensusScores_SummaryTable.txt", emit : consensus_summary
		path "ConsensusScores_SummaryTableByGroup.txt", emit : consensus_group_summary
		path "QC_composition.pdf", emit : composition_pdf
		path "CNV-Quadrants.pdf", emit : cnv_quadrants_pdf
		path "QC_composition_mqc.jpg", emit: composition_jpg
		path "CNV-Quadrants_mqc.jpg", emit: cnv_quadrants_jpg

	script:
	"""
	Rscript /usr/local/bin/cnvSummarizer.R \
		--rds_file ${rds_file} \
		--out_file AllSample-GinkgoSegmentSummary.txt
	echo "Finished cnvSummarizer.R"

	head AllSample-GinkgoSegmentSummary.txt
	echo "================================"

	# Check if selected_metrics is the dummy empty file
	if [ "${selected_metrics}" != "empty_selected_metrics.txt" ]; then
		SELECTED_METRICS_ARG="--selected_metrics ${selected_metrics}"
	else
		SELECTED_METRICS_ARG=""
	fi

	Rscript /usr/local/bin/function_plot_qc_dna.R \
		--seg_copy_file ${seg_copy_file} \
		--metrics_file ${metrics_file} \
		\$SELECTED_METRICS_ARG \
		--metadata_file ${metadata_file} \
		--cnv_summary_file AllSample-GinkgoSegmentSummary.txt \
		--plot_qc_config ${plot_qc_config}
	echo "Finished function_plot_qc_dna.R"

	Rscript /usr/local/bin/function_cnv_quadrants_qc.R \
		--metrics_file nf-preseq-pipeline_all_metrics_mqc_withcnv.txt
	echo "Finished function_cnv_quadrants_qc.R"

	mv nf-preseq-pipeline_all_metrics_mqc_withcnv.txt nf-preseq-pipeline_all_metrics_mqc.txt
	
	if [ -f "nf-preseq-pipeline_selected_metrics_mqc_withcnv.txt" ]; then
		mv nf-preseq-pipeline_selected_metrics_mqc_withcnv.txt nf-preseq-pipeline_selected_metrics_mqc.txt
	fi
	"""
}

process RNA_QC_PLOTS {   
    tag "qc_plots"
	publishDir "${publish_dir}_${params.timestamp}/tertiary_analyses/qc_plots", enabled:"$enable_publish"

	input:
		path matrix_file
		path metrics_file
		path metadata_file
		// path plot_qc_config
		val(publish_dir)
		val(enable_publish)

	output:
		path "composition_rnaqc.pdf"
		path "summary_verdict.txt"
		path "summary_verdict_group.txt"
		path "RNA-QC_ConsensusScores.txt", emit : consensus_scores
		path "RNAQC_composition_mqc.jpg", emit: composition_jpg

	script:
	"""
	Rscript /usr/local/bin/rna_qc_plot.R \
		--matrix_file ${matrix_file} \
		--metrics_file ${metrics_file} \
		--metadata_file ${metadata_file}
	echo "Finished rna_qc_plot script"
	"""
}
workflow QC_PLOTS_WF {
	take:
		ch_rds_file
		ch_seg_copy_file
		ch_metrics_file
		ch_selected_metrics
		ch_metadata_file
		ch_config_file
		ch_publish_dir
		ch_enable_publish

	main:
		QC_PLOTS(
			ch_rds_file,
			ch_seg_copy_file,
			ch_metrics_file,
			ch_selected_metrics,
			ch_metadata_file,
			ch_config_file,
			ch_publish_dir,
			ch_enable_publish
		)

	emit:
		ginkgo_summary = QC_PLOTS.out[0]
		allmetrics_with_cnv = QC_PLOTS.out.allmetrics_with_cnv
		selectedmetrics_with_cnv = QC_PLOTS.out.selectedmetrics_with_cnv
		consensus_scores = QC_PLOTS.out.consensus_scores
		consensus_summary = QC_PLOTS.out.consensus_summary
		consensus_group_summary = QC_PLOTS.out.consensus_group_summary
		composition_jpg = QC_PLOTS.out.composition_jpg
		cnv_quadrants_jpg = QC_PLOTS.out.cnv_quadrants_jpg
}

workflow RNA_QC_PLOTS_WF {
	take:
		ch_matrix_file
		ch_metrics_file
		ch_metadata_file
		ch_publish_dir
		ch_enable_publish

	main:
		RNA_QC_PLOTS(
			ch_matrix_file,
			ch_metrics_file,
			ch_metadata_file,
			ch_publish_dir,
			ch_enable_publish
		)

	emit:
		composition_pdf = RNA_QC_PLOTS.out[0]
		summary_verdict = RNA_QC_PLOTS.out[1]
		summary_verdict_group = RNA_QC_PLOTS.out[2]
		composition_jpg = RNA_QC_PLOTS.out.composition_jpg
}

workflow {
	if (params.qc_type != "dna") {
		RNA_QC_PLOTS_WF(
			params.matrix_file,
			params.rna_metrics_file,
			params.rna_metadata_file,
			params.publish_dir,
			params.enable_publish
		)
	} else {
		// Default to DNA QC plots
		// Create an empty file for selected_metrics when not provided
		def empty_selected_metrics = file('empty_selected_metrics.txt')
		empty_selected_metrics.text = ''
		
		QC_PLOTS_WF(
			params.rds_file,
			params.seg_copy_file,
			params.metrics_file,
			empty_selected_metrics,
			params.metadata_file,
			params.plot_qc_config,
			params.publish_dir,
			params.enable_publish
		)
	}
}
