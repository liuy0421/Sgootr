configfile: "config.yaml"
ruleorder: iteration_setup > prune

PATIENTS=list(config['PATIENTS'].keys())
OUTDIR = config['OUTDIR']

rule all:
    input:
        *[expand(OUTDIR+"{patient}/t{i}/{labels}.png",  patient=PATIENT, 
                                                        i=range(config['MAX_ITER']+1),
                                                        labels=list(config['PATIENTS'][PATIENT]['palette'].keys()))
          for PATIENT in PATIENTS]


# CNA input is optional
def get_data_setup_input(ws):
    if config['PATIENTS'][ws.patient]['cna'] == None:
        return config['PATIENTS'][ws.patient]['methylated_rc'], \
               config['PATIENTS'][ws.patient]['unmethylated_rc']
    return config['PATIENTS'][ws.patient]['methylated_rc'], \
           config['PATIENTS'][ws.patient]['unmethylated_rc'], \
           config['PATIENTS'][ws.patient]['cna']

rule data_setup:
    input:
        get_data_setup_input
    output:
        temp(OUTDIR+"{patient}/rc_cna.npz")
    params:
        default_cna = config['DEFAULT_CNA']
    log:
        OUTDIR+"{patient}/logs/data_setup.log"
    script:
        "scripts/data_setup.py"


rule error_correction:
    input:
        OUTDIR+"{patient}/rc_cna.npz"
    output:
        temp(OUTDIR+"{patient}/corrected_rc_cna.npz")
    params:
        e = config['E']
    threads: 
        4
    log:
        OUTDIR+"{patient}/logs/error_correction.log"
    script:
        "scripts/error_correction.py"


def get_post_filter_input(ws):

    return OUTDIR+ws.patient+"/corrected_rc_cna.npz", \
           config['PATIENTS'][ws.patient]['whitelist']

rule post_filter:
    input:
        get_post_filter_input
    output:
        OUTDIR+"{patient}/input.npz"
    params: 
        site_coverage_threshold = config['SITE_COVERAGE_THRESHOLD'],
        cell_coverage_threshold = config['CELL_COVERAGE_THRESHOLD']
    log:
        OUTDIR+"{patient}/logs/post_filter.log"
    script:
        "scripts/post_filter.py"


rule iteration_setup:
    input: 
        OUTDIR+"{patient}/input.npz"
    output:
        OUTDIR+"{patient}/t0/site_mask.npz",
        OUTDIR+"{patient}/heuristically_called_statuses.npz",
        OUTDIR+"{patient}/status_likelihoods.npz"
    threads:
        4
    params:
        status_confidence_threshold = config['STATUS_CONFIDENCE_THRESHOLD'],
        p00 = config['P00'],
        p10 = config['P10'],
        p11 = config['P11'],
        p = config['P']
    log:
        OUTDIR+"{patient}/logs/iteration_setup.log"
    script:
        "scripts/iteration_setup.py"


rule compute_distances:
    input:
        OUTDIR+"{patient}/t{i}/site_mask.npz",
        OUTDIR+"{patient}/status_likelihoods.npz"
    output:
        OUTDIR+"{patient}/t{i}/pairwise_distances.npz"
    params:
        d0010 = config['D0010'],
        d0011 = config['D0011'],
        d1011 = config['D1011'],
    threads:
        4
    log:
        OUTDIR+"{patient}/logs/t{i}/compute_distances.log"
    script:
        "scripts/compute_distances.py"


def get_build_tree_input(ws):

    i = int(ws.i)
    nwks = [OUTDIR+ws.patient+"/t{}/tree.nwk".format(round) for round in range(i)]    

    if i == 0:
        return OUTDIR+ws.patient+"/t0/pairwise_distances.npz"
    elif i > 0:
        return OUTDIR+ws.patient+"/t{}/pairwise_distances.npz".format(i), \
               OUTDIR+ws.patient+"/t{}/RF.txt".format(i-1), \
               *nwks

rule build_tree:
    input:
        get_build_tree_input
    output:
        OUTDIR+"{patient}/t{i}/tree.nwk",
        OUTDIR+"{patient}/t{i}/RF.txt"
    params:
        root = lambda ws: config['PATIENTS'][ws.patient]['root']
    log:
        OUTDIR+"{patient}/logs/t{i}/build_tree.log"
    conda:
        "envs/skbio.yaml"
    shell:
        "python scripts/build_tree.py {output} {params.root} {log} {input}"


def get_prune_input(ws):

    i = int(ws.i)
    nwks = [OUTDIR+ws.patient+"/t{}/tree.nwk".format(round) for round in range(i)]    

    return OUTDIR+ws.patient+"/t{}/tree.nwk".format(i-1), \
           OUTDIR+ws.patient+"/t{}/site_mask.npz".format(i-1), \
           OUTDIR+ws.patient+"/heuristically_called_statuses.npz"

rule prune:
    input:
        get_prune_input
    output:
        OUTDIR+"{patient}/t{i}/site_mask.npz",
        OUTDIR+"{patient}/t{i}/persistence_scores.npz"
    params:
        kappa = config['KAPPA'],
        partition_validity_threshold = config['PARTITION_VALIDITY_THRESHOLD'],
        minimum_subtree_size = config['MINIMUM_SUBTREE_SIZE']
    threads:
        4
    log:
        OUTDIR+"{patient}/logs/t{i}/prune.log"
    conda:
        "envs/skbio.yaml"
    shell:
        "python scripts/prune.py {input} {output} {params.kappa} " \ 
        "{params.partition_validity_threshold} {params.minimum_subtree_size} " \
        "{threads} {log}"


def get_visualize_input(ws):

    return OUTDIR+"{}/t{}/tree.nwk".format(ws.patient, ws.i), \
           config['PATIENTS'][ws.patient]["labels"]

rule visualize:
    input:
        get_visualize_input
    output:
        OUTDIR+"{patient}/t{i}/{labels}.png"        
    params:
        palette = lambda ws: config['PATIENTS'][ws.patient]['palette'][ws.labels],
        color_by = "{labels}" 
    script:
        "scripts/visualize_tree.R"



