configfile: "config.yaml"
ruleorder: iteration_setup > prune

# CNA input is optional
def get_data_setup_input(ws):
    if config[ws.patient]['cna'] == None:
        return config[ws.patient]['methylated_rc'], \
               config[ws.patient]['unmethylated_rc']
    return config[ws.patient]['methylated_rc'], \
           config[ws.patient]['unmethylated_rc'], \
           config[ws.patient]['cna']

rule data_setup:
    input:
        get_data_setup_input
    output:
        temp("{patient}/rc_cna.npz")
    params:
        default_cna = config['DEFAULT_CNA']
    log:
        "{patient}/logs/data_setup.log"
    script:
        "scripts/data_setup.py"


rule error_correction:
    input:
        "{patient}/rc_cna.npz"
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
        "{patient}/t0/site_mask.npz",
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


rule compute_distances:
    input:
        "{patient}/t{i}/site_mask.npz",
        "{patient}/status_likelihoods.npz"
    output:
        "{patient}/t{i}/pairwise_distances.npz"
    params:
        d0010 = config['D0010'],
        d0011 = config['D0011'],
        d1011 = config['D1011'],
    threads:
        4
    log:
        "{patient}/logs/t{i}/compute_distances.log"
    script:
        "scripts/compute_distances.py"


def get_build_tree_input(ws):

    i = int(ws.i)
    nwks = [ws.patient+"/t{}/tree.nwk".format(round) for round in range(i)]    

    if i == 0:
        return ws.patient+"/t0/pairwise_distances.npz"
    elif i > 0:
        return ws.patient+"/t{}/pairwise_distances.npz".format(i), \
               ws.patient+"/t{}/RF.txt".format(i-1), \
               *nwks

rule build_tree:
    input:
        get_build_tree_input
    output:
        "{patient}/t{i}/tree.nwk",
        "{patient}/t{i}/RF.txt"
    params:
        root = lambda ws: config[ws.patient]['root']
    log:
        "{patient}/logs/t{i}/build_tree.log"
    conda:
        "gmelin-larch"
    shell:
        "python scripts/build_tree.py {output} {params.root} {log} {input}"


def get_prune_input(ws):

    i = int(ws.i)
    nwks = [ws.patient+"/t{}/tree.nwk".format(round) for round in range(i)]    

    return ws.patient+"/t{}/tree.nwk".format(i-1), \
           ws.patient+"/t{}/site_mask.npz".format(i-1), \
           ws.patient+"/heuristically_called_statuses.npz"

rule prune:
    input:
        get_prune_input
    output:
        "{patient}/t{i}/site_mask.npz"
    params:
        kappa = config['KAPPA'],
        partition_validity_threshold = config['PARTITION_VALIDITY_THRESHOLD']
    threads:
        4
    log:
        "{patient}/logs/t{i}/prune.log"
    conda:
        "gmelin-larch"
    shell:
        "python scripts/prune.py {input} {output} {params.kappa} " \ 
        "{params.partition_validity_threshold} {threads} {log}"


def get_visualize_input(ws):

    return "{}/t{}/tree.nwk".format(ws.patient, ws.i), \
           config[ws.patient]["labels"]

rule visualize:
    input:
        get_visualize_input
    output:
        "{patient}/t{i}/{labels}.png"        
    params:
        palette = lambda ws: config[ws.patient]['palette'][ws.labels],
        color_by = "{labels}" 
    script:
        "scripts/visualize_tree.R"



