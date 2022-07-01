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
def get_qc_input(ws):
    if config[ws.patient]['cna'] == None:
        return config[ws.patient]['methylated_rc'], \
               config[ws.patient]['unmethylated_rc']
    return config[ws.patient]['methylated_rc'], \
           config[ws.patient]['unmethylated_rc'], \
           config[ws.patient]['cna']

rule qc:
    input:
        get_qc_input
    output:
        temp("{patient}/qc_rc_cna.npz")
    params:
        default_cna = config['DEFAULT_CNA'],
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
        e = config['E']
    threads: 
        4
    log:
        "{patient}/logs/error_correction.log"
    script:
        "scripts/error_correction.py"


def get_post_filter_input(ws):

    return "{}/corrected_rc_cna.npz".format(ws.patient), \
           config[ws.patient]['whitelist']

rule post_filter:
    input:
        get_post_filter_input
    output:
        "{patient}/input.npz"
    params: 
        site_coverage_threshold = config['SITE_COVERAGE_THRESHOLD'],
        cell_coverage_threshold = config['CELL_COVERAGE_THRESHOLD']
    log:
        "{patient}/logs/post_filter.log"
    script:
        "scripts/post_filter.py"


rule iteration_setup:
    input: 
        "{patient}/input.npz"
    output:
        "{patient}/t0_site_mask.npz",
        "{patient}/heuristically_called_statuses.npz",
        "{patient}/status_likelihoods.npz"
    threads:
        4
    params:
        status_confidence_threshold = config['STATUS_CONFIDENCE_THRESHOLD'],
        p00 = config['P00'],
        p10 = config['P10'],
        p11 = config['P11'],
        p = config['P']
    log:
        "{patient}/logs/iteration_setup.log"
    script:
        "scripts/iteration_setup.py"


def get_compute_distances_input(ws):

    i = int(ws.i)

    if i == 0:
        return ws.patient+"/heuristically_called_statuses.npz", \
               ws.patient+"/status_likelihoods.npz", \
               ws.patient+"/site_mask.npz"
    elif i > 0:
        return ws.patient+"/heuristically_called_statuses.npz", \
               ws.patient+"/status_likelihoods.npz", \
               ws.patient+"/site_mask.npz", \
               ws.patient+"/persistence_scores_t{}.npz".format(i-1), \
               ws.patient+"/t{}.nwk".format(i-1), \
               ws.patient+"/RF_{}.txt".format(i-1)
    else:
        raise ValueError('Invalid value for number of iterations')


rule compute_distances:
    input:
        get_compute_distances_input
    output:
        "{patient}/t{i}_pairwise_distances.npz"#,
    params:
        d0010 = config['D0010'],
        d0011 = config['D0011'],
        d1011 = config['D1011'],
    threads:
        4
    log:
        "{patient}/logs/compute_distances_{i}.log"
    script:
        "scripts/compute_distances.py"


rule build_tree:
    input:
        "{patient}/t{i}_pairwise_distances.npz"
    output:
        "{patient}/t{i}.nwk"
    params:
        root = lambda ws: config[ws.patient]['root']
    log:
        "{patient}/logs/build_tree_{i}.log"
    conda:
        #"envs/skbio.yaml"
        "gmelin-larch"
    #script:
    #    "scripts/script_build_tree.py"
    shell:
        "python scripts/build_tree.py {input} {output} {params.root} {log}"



