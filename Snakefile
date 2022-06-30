configfile: "config.yaml"
'''
rule all:
    input:
        dynamic(expand("{patient}/iterations/scores{i}.npz", patient=config.keys)),
        dynamic(expand("{patient}/iterations/t{i}.nwk", patient=config.keys)),
        expand("{patient}/final.nwk", patient=config.keys),
        expand("{patient}/sites.npz", patient=config.keys)
'''
# CNA input is optional
def get_qc_input(wildcards):
    if config[wildcards.patient]['cna'] == None:
        return config[wildcards.patient]['methylated_rc'], \
               config[wildcards.patient]['unmethylated_rc']
    return config[wildcards.patient]['methylated_rc'], \
           config[wildcards.patient]['unmethylated_rc'], \
           config[wildcards.patient]['cna']

rule qc:
    input:
        get_qc_input
    output:
        temp("{patient}/qc_rc_cna.npz")
    params:
        default_cna=2,
        coverage_threshold=.66
    log:
        "{patient}/logs/qc.log"
    script:
        "scripts/qc.py"


rule error_correction:
    input:
        "{patient}/qc_rc_cna.npz"
    output:
        "{patient}/corrected_rc_cna.npz"
    params:
        e=.01
    threads: 
        4
    log:
        "{patient}/logs/error_correction.log"
    script:
        "scripts/error_correction.py"


def get_post_filter_input(wildcards):

    return "{}/corrected_rc_cna.npz".format(wildcards.patient), \
           config[wildcards.patient]['whitelist']

rule post_filter:
    input:
        get_post_filter_input
    output:
        "{patient}/input.npz"
    params: 
        site_coverage_threshold=.66,
        cell_coverage_threshold=.5
    log:
        "{patient}/logs/post_filter.log"
    script:
        "scripts/post_filter.py"


rule iteration_setup:
    input: 
        "{patient}/input.npz"
    output:
        "{patient}/site_mask.npz",
        "{patient}/heuristically_called_statuses.npz",
        "{patient}/status_likelihoods.npz"
    threads:
        4
    params:
        status_confidence_threshold=2,
        p00=.33,
        p10=.33,
        p11=.33,
        p=.5
    log:
        "{patient}/logs/iteration_setup.log"
    script:
        "scripts/iteration_setup.py"

'''
rule prune:
    input:
        "{patient}/input.npz"
    output:
        dynamic("{patient}/iterations/scores{i}.npz"),
        dynamic("{patient}/iterations/t{i}.nwk"),
        "{patient}/final.nwk",
        "{patient}/sites.npz"
    params:
        kappa=.1,
        max_iter=20,
        root=NC_507
    threads:
        4
    log:
        "{patient}/logs/prune.log"
    script:
        "scripts/prune.py"

'''
