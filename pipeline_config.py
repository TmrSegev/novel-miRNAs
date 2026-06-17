"""Shared configuration for per-library miRNA discovery pipelines."""

import os

DEFAULT_BASE_PATH = "/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq"
DEFAULT_NCRNA_DIR = "/sise/vaksler-group/IsanaRNA/Isana_Tzah/RNAcentral/ncRNAs_Caenorhabditis"

SPECIES_CONFIG = {
    "Elegans": {
        "libraries": [
            "CE57", "CE58", "CE59", "CE60", "CE61", "CE62",
            "CE63", "CE69", "CE78", "CE79", "CE80", "CE81",
        ],
        "srnabench_out_subdir": "sRNAtoolboxDB/out/Elegans",
        "mirdeep_out_subdir": "Elegans/mirdeep_out",
        "scripts_subdir": "Elegans/scripts",
        "good_candidates_subdir": "Elegans/good_candidates",
        "srnabench_folder_prefix": "Elegans_",
        "support_mode": "distinct_libraries",
    },
    "Macrosperma": {
        "libraries": ["MR4", "MR5", "MR6", "MR7", "MR8"],
        "srnabench_out_subdir": "sRNAtoolboxDB/out/Macrosperma",
        "mirdeep_out_subdir": "Macrosperma/mirdeep_out",
        "scripts_subdir": "Macrosperma/scripts",
        "good_candidates_subdir": "Macrosperma/good_candidates",
        "srnabench_folder_prefix": "Macrosperma_",
        "support_mode": "distinct_libraries",
    },
    "Sulstoni": {
        "libraries": ["SR0", "SR1", "SR2", "SR3", "SR4", "SR5", "SR6", "SR7"],
        "srnabench_out_subdir": "sRNAtoolboxDB/out/Sulstoni",
        "mirdeep_out_subdir": "Sulstoni/mirdeep_out",
        "scripts_subdir": "Sulstoni/scripts",
        "good_candidates_subdir": "Sulstoni/good_candidates",
        "srnabench_folder_prefix": "Sulstoni_",
        "support_mode": "distinct_libraries",
    },
    "Hofstenia": {
        "libraries": [
            "EC1", "EC2", "EC3", "GA1", "GA2", "GA3", "DI1", "DI2", "DI3",
            "PDi1", "PDi2", "PDi3", "PDii1", "PDii2", "PDii3", "PL1", "PL2", "PL3",
            "PH1", "PH2", "PH3", "HL1", "HL2", "HL3", "IST1", "IST2", "IST3",
            "AMP1", "AMP2", "AMP3", "SMA1", "SMA2", "SMA3",
        ],
        "srnabench_out_subdir": "sRNAtoolboxDB/out",
        "mirdeep_out_subdir": "Hofstenia/mirdeep_out",
        "scripts_subdir": "Hofstenia/scripts",
        "good_candidates_subdir": "Hofstenia/good_candidates",
        "srnabench_folder_prefix": "Hofstenia_",
        "support_mode": "dev_condition_replicates",
    },
}


def get_species_config(species, base_path=None):
    if species not in SPECIES_CONFIG:
        raise ValueError(
            f"Unknown species '{species}'. "
            f"Choose from: {', '.join(SPECIES_CONFIG.keys())}"
        )
    cfg = SPECIES_CONFIG[species].copy()
    cfg["species"] = species
    cfg["base_path"] = base_path or DEFAULT_BASE_PATH
    cfg["scripts_dir"] = os.path.join(cfg["base_path"], cfg["scripts_subdir"])
    cfg["good_candidates_dir"] = os.path.join(cfg["base_path"], cfg["good_candidates_subdir"])
    cfg["srnabench_out_dir"] = os.path.join(cfg["base_path"], cfg["srnabench_out_subdir"])
    cfg["mirdeep_out_dir"] = os.path.join(cfg["base_path"], cfg["mirdeep_out_subdir"])
    return cfg


def srnabench_folder_name(cfg, library):
    return f"{cfg['srnabench_folder_prefix']}{library}"


def mirdeep_folder_name(cfg, library):
    return library


def dev_condition(cfg, library_name, tool_name):
    """Return grouping key for good_candidates support filtering."""
    lib = library_name.split("_")[-1]
    if cfg["support_mode"] == "dev_condition_replicates":
        if tool_name == "sRNAbench":
            return lib[:-1] if len(lib) > 1 else lib
        return lib[:-1] if len(lib) > 1 else lib
    return lib


def build_description(series_cols, include_mirbase=False):
    """Unified Description: __ join, ; -> |, strip ID= and periods."""
    if include_mirbase:
        cols = ["Description_mirdeep", "Description_sRNAbench", "Description_mirbase"]
    else:
        cols = ["Description_mirdeep", "Description_sRNAbench"]
    return (
        series_cols[cols]
        .astype(str)
        .agg("__".join, axis=1)
        .str.replace(";", "|", regex=False)
        .str.replace("ID=", "", regex=False)
        .str.replace(".", "", regex=False)
        .str.replace("nan__", "", regex=False)
        .str.replace("__nan", "", regex=False)
    )
