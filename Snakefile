import re


## TODO: Implement check to see if conda environment has been created and activated

configfile: 'config.yaml'

sample_names = [re.sub('\..*', '', sample) for sample in config['samples']]
nthread = config['nthread']
reference = config['reference']
segmenter = config['segmenter']


## TODO: Add segmenter name to output segmentation files

# ---- 0. Get cnvkit.py scripts from GitHub and make symlink in top level directory
rule pull_cnvkit_github:
    output:
        'bin/cnvkit/cnvkit.py'
    shell:
        'git clone https://github.com/etal/cnvkit bin/cnvkit;'
        'ln -s bin/cnvkit/cnvkit.py .'


# ---- 1. Preprocess target and antitarget bins
rule equalize_target_bins:
    input:
        expand('metadata/{targets}.bed', targets=config['targets'])
    output:
        'metadata/binnedTargets.bed'
    shell:
        'python3 cnvkit.py target {input} --split -o {output}'


rule autocalculate_offtarget_bins:
    input:
        samples=expand('procdata/{samples}.bam', samples=config['samples']),
        target='metadata/binnedTargets.bed'
    output:
        'results/binnedTargets.target.bed',
        'results/binnedTargets.antitarget.bed'
    shell:
        'python3 cnvkit.py autobin {input.samples} -t {input.target} --method amplicon;'
        'mv *.bed results;'


rule calculate_per_sample_target_antitarget_coverage:
    input:
        samples=expand('procdata/{samples}.bam', samples=config['samples']),
        target_bed='results/binnedTargets.target.bed',
        antitarget_bed='results/binnedTargets.antitarget.bed'
    output:
        target=expand('results/{sample_names}.targetcoverage.cnn', sample_names=sample_names),
        antitarget=expand('results/{sample_names}.antitargetcoverage.cnn', sample_names=sample_names)
    run:
        for i in range(len(input.samples)):
            sample, target, antitarget = input.samples[i], output.target[i], output.antitarget[i]
            shell(f'python3 cnvkit.py coverage {sample} {input.target_bed} -o {target} -p {nthread}')
            shell(f'python3 cnvkit.py coverage {sample} {input.antitarget_bed} -o {antitarget} -p {nthread}')


# ---- 2. Create copy number calls for reference
## If this step breaks make sure to delete the reference fai file, it could be corrupted
rule call_reference_copy_numbers:
     input:
        reference=expand('reference/{reference}.fa', reference=reference),
        targets='results/binnedTargets.target.bed'
     output:
        expand('reference/{reference}.cnn', reference=reference)
     shell:
        'python3 cnvkit.py reference -o {output} -f {input.reference} -t {input.targets}'


# ---- 3. Calculate per sample copy number ratio, adjusting for regional coverage and GC bias
rule call_sample_copynumber_ratios:
    input:
        sample_targets=expand('results/{sample_names}.targetcoverage.cnn', sample_names=sample_names),
        sample_antitargets=expand('results/{sample_names}.antitargetcoverage.cnn', sample_names=sample_names),
        reference=expand('reference/{reference}.cnn', reference=reference)
    output:
        results=expand('results/{sample_names}.cnr', sample_names=sample_names)
    run:
        for i in range(len(input.sample_targets)):
            target, antitarget, result = input.sample_targets[i], input.sample_antitargets[i], output.results[i]
            shell('python3 cnvkit.py fix {target} {antitarget} {input.reference} -o {result}')


# ---- 4. Segment per sample copy number ratios, dropping outliers and low coverage regions
rule segment_sample_copynumber_ratios:
    input:
        cn_ratios=expand('results/{sample_names}.cnr', sample_names=sample_names)
    output:
        results=expand('results/{sample_names}.cns', sample_names=sample_names)
    run:
        for i in range(len(input.cn_ratios)):
            ratio, result = input.cn_ratios[i], output.results[i]
            shell('python3 cnvkit.py segment {ratio} -o {result} -m {segmenter} --drop-low-coverage --drop-outliers')


# ---- 5. Call integer copy numbers for each segment
rule call_sample_integer_copynumbers:
    input:
        expand('results/{sample_names}.cns', sample_names=sample_names)
    output:
        expand('results/{sample_names}.call.cns', sample_names=sample_names) # Including output for snakemake file checking
    run:
        for i in range(len(input)):
            segment, result = input[i], output[i]
            shell('python3 cnvkit.py call {segment} -y -m threshold -t=-1.1,-0.4,0.3,0.7 -o {result}')


# ---- 6. Generate per sample copy number ratio/segment scatter plots
rule scatterplot_sample_copynumber:
    input:
        cns=expand('results/{sample_names}.cns', sample_names=sample_names),
        cnr=expand('results/{sample_names}.cnr', sample_names=sample_names)
    output:
        expand('results/figures/{sample_names}_scatter.svg', sample_names=sample_names)
    run:
        for i in range(len(input.cns)):
            cns, cnr, result = input.cns[i], input.cnr[i], output[i]
            shell('python3 cnvkit.py scatter {cnr} -s {cns} -o {result}')


# ---- 7. Generate all sample heatmap
rule heatmap_copynumber_summary:
    input:
        expand('results/{sample_names}.cns', sample_names=sample_names)
    output:
        'results/figures/lmsHaloplexHeatmap.pdf'
    shell:
        'python3 cnvkit.py heatmap {input} -o {output}'


# ---- 8. Generate per sample chrosome diagrams
rule plot_sample_chromosome_diagrams:
    input:
        cns=expand('results/{sample_names}.cns', sample_names=sample_names),
        cnr=expand('results/{sample_names}.cnr', sample_names=sample_names)
    output:
        expand('results/figures/{sample_names}_chromosome.pdf', sample_names=sample_names)
    threads: 14
    run:
        for i in range(len(input.cnr)):
            cns, cnr, result = input.cns[i], input.cnr[i], output[i]
            shell('python3 cnvkit.py diagram {cnr} -s {cns} -o {result}')

#https://harvardlms.blob.core.windows.net/lms/LMSARC/harvardlms/LMSmodelSystems/CancerCellLines/HaloplexAnalysis/