print(
    """
    \tAlignment Pipeline 2.0 (Paired-end FASTQ / BAM / SRA)
    \tAuthor: Shashwat Sahay
    \tEmail: shashwat.sahay@charite.de
    \tVersion: 2.0.0

    """
)


def get_data_time():
    now = datetime.now()

    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    return "========= " + dt_string + " ========= "


print(get_data_time(), "Starting run")

#######################################################################
################# Sanity check on the input data ######################
#######################################################################

if "sample" not in config:
    raise ValueError("Sample name not provided")

if "pid" not in config:
    raise ValueError("Patient ID not provided")

_metadata_cfg = config.get("metadata")
_bamfile_cfg = config.get("bamfile")
_sra_cfg = config.get("sra") or config.get("sra_id")
ngc_file = config.get("ngc") or config.get("ngc_file") or None

_provided_inputs = [
    name
    for name, val in [("metadata", _metadata_cfg), ("bamfile", _bamfile_cfg), ("sra", _sra_cfg)]
    if val
]

if len(_provided_inputs) == 0:
    raise ValueError("Provide one of 'metadata', 'bamfile', or 'sra' in config.")
if len(_provided_inputs) > 1:
    raise ValueError(
        f"Provide only one of 'metadata', 'bamfile', or 'sra' in config (got: {', '.join(_provided_inputs)})."
    )

if _metadata_cfg:
    input_mode = "fastq"
    bamfile = None
    sra_id = None
    sep = ","
    if _metadata_cfg.split(".")[-1] in ["tsv", "txt"]:
        sep = "\t"
    elif _metadata_cfg.split(".")[-1] not in ["csv", "tsv", "txt"]:
        raise ValueError(
            "Can't resolve the separator for the metadata file from extension"
        )
    metadata = pd.read_csv(_metadata_cfg, dtype=str, sep=sep)
    metadata = metadata[
        (metadata["SAMPLE_TYPE"] == config["sample"])
        & (metadata["PATIENT_ID"] == config["pid"])
    ]
    if metadata.shape[0] == 0:
        raise ValueError(
            "No metadata found for the given sample: %s and patient ID: %s"
            % (config["sample"], config["pid"])
        )
    print(get_data_time(), f"Input mode: FASTQ metadata ({_metadata_cfg})")
elif _bamfile_cfg:
    input_mode = "bam"
    metadata = None
    bamfile = _bamfile_cfg
    sra_id = None
    print(get_data_time(), f"Input mode: BAM ({bamfile})")
else:
    input_mode = "sra"
    metadata = None
    bamfile = None
    sra_id = str(_sra_cfg).strip()
    print(get_data_time(), f"Input mode: SRA ({sra_id})")
    if ngc_file:
        if not os.path.exists(ngc_file):
            raise ValueError(f"Provided ngc file does not exist: {ngc_file}")
        print(get_data_time(), f"Using NGC key file: {ngc_file}")

##########################################################
###### Setting Sequencing type and related parameters ####
##########################################################

if "SeqType" not in config:
    raise ValueError("Sequencing type not provided")
else:
    seq_type = config["SeqType"]
    if seq_type not in ["WGS", "WES", "Panel"]:
        raise ValueError("Invalid sequencing type provided")
    if seq_type in ["WES", "Panel"]:
        if "target_regions" not in config:
            raise ValueError("Target regions not provided")
        else:
            target_regions = config["target_regions"]

        if "chrom_sizes" not in config:
            raise ValueError("Chromosome sizes not provided")
        else:
            chrom_sizes = config["chrom_sizes"]
        if "dict_genome" not in config:
            raise ValueError("Genome not provided")
        else:
            dict_genome = config["dict_genome"]

        if "bait_regions" not in config:
            print(
                get_data_time(),
                "Bait regions not provided will be calculated from target regions",
            )
        else:
            bait_regions = config["bait_regions"]


def _parse_adapter_list(val):
    if isinstance(val, list):
        return val
    return [x.strip() for x in str(val).replace(";", ",").split(",") if x.strip()]


def resolve_sample_param(config_key, meta_cols, default_val, as_list=False, required=False):
    if metadata is not None:
        for col in meta_cols:
            if col in metadata.columns:
                vals = [
                    str(v).strip()
                    for v in metadata[col].dropna().tolist()
                    if str(v).strip() != ""
                ]
                unique_vals = list(dict.fromkeys(vals))
                if len(unique_vals) > 1:
                    raise ValueError(
                        f"Multiple conflicting values for '{col}' in metadata for sample: {unique_vals}"
                    )
                if len(unique_vals) == 1:
                    return _parse_adapter_list(unique_vals[0]) if as_list else unique_vals[0]

    cfg_val = config.get(config_key)
    if cfg_val is not None and cfg_val != "" and cfg_val != []:
        return _parse_adapter_list(cfg_val) if as_list else cfg_val

    if required:
        raise ValueError(f"'{config_key}' not provided in metadata or config")

    print(
        get_data_time(),
        f"'{config_key}' not provided in metadata or config; defaulting to {repr(default_val)}",
    )
    return default_val


##########################################################
################ Setting Library kit #####################
##########################################################

library_prep_kit = resolve_sample_param(
    "library_prep_kit", ["LIBRARY_PREP_KIT", "library_prep_kit"], "Unknown"
)


##########################################################
################## Setting Genome ########################
##########################################################

if "genome" not in config:
    raise ValueError("Genome not provided")
else:
    genome = config["genome"]

if "dbsnp" not in config:
    raise ValueError("dbSNP file not provided")
else:
    dbsnp = config["dbsnp"]


##########################################################
############# Setting working directory ##################
##########################################################

if "work_dir" not in config:
    raise ValueError("Work directory not provided")
else:
    wrkdir = Path(config["work_dir"])

if "log_dir" not in config or not config["log_dir"]:
    logdir = wrkdir / "logs"
else:
    logdir = Path(config["log_dir"])

if "scratch_dir" not in config or not config["scratch_dir"]:
    scratch_path = wrkdir
else:
    scratch_path = Path(config["scratch_dir"])
    if not os.path.exists(scratch_path):
        print(get_data_time(), "Scratch directory does not exist, creating it")
        os.makedirs(scratch_path, exist_ok=True)

scratch_dir = str(scratch_path)

print(get_data_time(), "Setting working directory to %s" % wrkdir)
print(get_data_time(), "Setting log directory to %s" % logdir)
print(get_data_time(), "Setting temp directory to %s" % scratch_dir)


#########################################################################
############# Setting variables related to trimming #####################
#########################################################################

if "trim_adapters" not in config or (not config["trim_adapters"]):
    print(get_data_time(), "Adapter trimming has been turned off")
    config["trim_adapters"] = False
    adapter_seq_r1 = None
    adapter_seq_r2 = None
else:
    adapter_seq_r1 = resolve_sample_param(
        "adapter_forward_read",
        ["ADAPTER_FORWARD_READ", "adapter_forward_read"],
        None,
        as_list=True,
        required=True,
    )
    adapter_seq_r2 = resolve_sample_param(
        "adapter_reverse_read",
        ["ADAPTER_REVERSE_READ", "adapter_reverse_read"],
        None,
        as_list=True,
        required=True,
    )


if "max_coverage" not in config:
    print(get_data_time(), "Setting default value for max_coverage to 10000")
    max_coverage = 10000
else:
    max_coverage = config["max_coverage"]


# Extract unique lanes and runs (from metadata in fastq mode, or single virtual run/lane for bam/sra)
if input_mode == "fastq":
    LANE = metadata["LANE_NO"].unique().tolist()
    RUN_ID = metadata["RUN_ID"].unique().tolist()
else:
    LANE = ["L0001"]
    RUN_ID = ["run1"]


def filter_combinator(combinator, allow_list):
    def filtered_combinator(*args, **kwargs):
        for wc_comb in combinator(*args, **kwargs):
            if frozenset(wc_comb) in allow_list:
                yield wc_comb

    return filtered_combinator


# Allowed combinations for sample, run, and lane
allow_list = set()
if input_mode == "fastq":
    for run_id in RUN_ID:
        for lane in metadata[metadata.RUN_ID == run_id]["LANE_NO"].unique().tolist():
            allow_list.add(
                frozenset(
                    {
                        ("run_id", run_id),
                        ("sample", config["pid"] + "_" + config["sample"]),
                        ("lane", lane),
                    }
                )
            )
else:
    allow_list.add(
        frozenset(
            {
                ("run_id", "run1"),
                ("sample", config["pid"] + "_" + config["sample"]),
                ("lane", "L0001"),
            }
        )
    )

filtered_product = filter_combinator(product, allow_list)

# Allowed combinations for FastQC on trimmed reads (only when trim_adapters is enabled)
allow_list_2 = set()
if config["trim_adapters"]:
    for read in ["R1", "R2"]:
        for comb in allow_list:
            allow_list_2.add(comb | frozenset({("read", read)}))

filtered_product_2 = filter_combinator(product, allow_list_2)

# ==============================================================================
# extra params
# ==============================================================================

extra_cutadapt_params = resolve_sample_param(
    "extra_cutadapt_params", ["EXTRA_CUTADAPT_PARAMS", "extra_cutadapt_params"], ""
)
extra_samtools_convert_params = resolve_sample_param(
    "extra_samtools_convert_params",
    ["EXTRA_SAMTOOLS_CONVERT_PARAMS", "extra_samtools_convert_params"],
    "",
)
mosdepth_target_regions = target_regions if seq_type in ["Panel", "WES"] else []
mosdepth_by_flag = f"--by {target_regions}" if seq_type in ["Panel", "WES"] else ""


def getFastq(wildcards, metadata, base_dir, trim_adapters):
    if trim_adapters:
        fastq_r1 = [
            base_dir
            / "fastq"
            / wildcards.run_id
            / "cutadapt"
            / (wildcards.sample + "_R1_" + wildcards.lane + "_trim.fastq.gz"),
            base_dir
            / "fastq"
            / wildcards.run_id
            / "cutadapt"
            / (wildcards.sample + "_R2_" + wildcards.lane + "_trim.fastq.gz"),
        ]
    elif input_mode == "fastq":
        fastq_r1 = (
            metadata[
                (metadata["RUN_ID"] == wildcards.run_id)
                & (metadata["LANE_NO"] == wildcards.lane)
                & (metadata["SAMPLE_NAME"] == wildcards.sample)
            ]
            .sort_values("READ")["FASTQ_FILE"]
            .values
        )
    else:
        fastq_r1 = [
            base_dir
            / "fastq"
            / "raw"
            / wildcards.run_id
            / (wildcards.sample + "_R1_" + wildcards.lane + ".fastq.gz"),
            base_dir
            / "fastq"
            / "raw"
            / wildcards.run_id
            / (wildcards.sample + "_R2_" + wildcards.lane + ".fastq.gz"),
        ]
    if len(fastq_r1) != 2:
        print(fastq_r1)
        raise ValueError("Expected two files R1 and R2")

    return fastq_r1
