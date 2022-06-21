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
        default_cna=2
    log:
        "{patient}/logs/qc.log"
    script:
        "scripts/qc.py"


rule error_correction:
    input:
        "{patient}/qc_rc_cna.npz"
    output:
        temp("{patient}/corrected_rc_cna.npz")
    params:
        e=.01
    threads: 
        4
    log:
        "{patient}/logs/error_correction.log"
    script:
        "scripts/error_correction.py"

'''
rule filter:
    input:
        "{patient}/corrected_rc_cna.npz"
    output:
        protected("{patient}/t0.npz")
    params: 
        homogeneity_pi=.98
    log:
        "{patient}/logs/filter.log"
    script:
        "scripts/filter.py"

rule prune:
    input:
        "{patient}/t0.npz"
    output:
        dynamic("{patient}/iterations/scores{i}.npz"),
        dynamic("{patient}/iterations/t{i}.nwk"),
        "{patient}/final.nwk",
        "{patient}/sites.npz"
    params:
        pruning_fraction=.05
    threads:
        64
    log:
        "{patient}/logs/prune.log"
    script:
        "scripts/prune.py"

'''
