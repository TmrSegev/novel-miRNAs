"""Shared configuration for per-library miRNA discovery pipelines."""

import os

DEFAULT_BASE_PATH = "/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq"
DEFAULT_NCRNA_DIR = "/mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/ncRNAs_Caenorhabditis"

# Global filter defaults (Hofstenia values; overridable via CLI).
# Historical nematode values: filter_mc=100, exclude_c=1000, low_score_total=1000
FILTER_MC_DEFAULT = 10
EXCLUDE_C_DEFAULT = 100
LOW_SCORE_TOTAL_DEFAULT = 100

SEED_FILE_NEMATODE = "mirbase_data/Seeds.txt"
SEED_FILE_HOFSTENIA = "mirbase_data/ALL_seed_family_from_mirgendb.csv"

DEFAULT_BLAST_DB = os.path.join(
    DEFAULT_BASE_PATH, "mirbase_data/Caenorhabditis_pre_miRNAsDB"
)

NEMATODE_PROFILE = {
    "seed_file": SEED_FILE_NEMATODE,
    "seed_sep": "\t",
    "seed_encoding": "utf-8",
    "use_blast": True,
    "use_mirbase_intersects": False,
    "ziv_profile": "structural_sheets",
    "mirge_input_sheet": "(D) Structural Features",
    "apply_manual_corrections": False,
    "compare_genome_qc": True,
    "run_sensitivity_plots": False,
    "ncrna_dir": DEFAULT_NCRNA_DIR,
    "output_subdir": "RNAcentral/miRNAs/{species}",
}

HOFSTENIA_PROFILE = {
    "seed_file": SEED_FILE_HOFSTENIA,
    "seed_sep": ",",
    "seed_encoding": "latin-1",
    "use_blast": False,
    "use_mirbase_intersects": False,
    "ziv_profile": "unfiltered_only",
    "mirge_input_sheet": "(A) Unfiltered",
    "apply_manual_corrections": True,
    "compare_genome_qc": False,
    "run_sensitivity_plots": False,
    "ncrna_dir": DEFAULT_NCRNA_DIR,
    "output_subdir": "RNAcentral/miRNAs/{species}",
}

HOFSTENIA_NEW_GENOME_VARIANT = {
    "mirdeep_out_subdir": "Hofstenia_newGenome/mirdeep_out",
    "scripts_subdir": "Hofstenia_newGenome/scripts",
    "good_candidates_subdir": "Hofstenia_newGenome/good_candidates",
    "srnabench_out_subdir": "sRNAtoolboxDB/out/Hofstenia_newGenome",
    "output_subdir": "RNAcentral/miRNAs/Hofstenia_newGenome",
    "species_label": "Hofstenia",
}

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
        **NEMATODE_PROFILE,
        "use_mirbase_intersects": True,
        "run_sensitivity_plots": True,
    },
    "Macrosperma": {
        "libraries": ["MR4", "MR5", "MR6", "MR7", "MR8"],
        "srnabench_out_subdir": "sRNAtoolboxDB/out/Macrosperma",
        "mirdeep_out_subdir": "Macrosperma/mirdeep_out",
        "scripts_subdir": "Macrosperma/scripts",
        "good_candidates_subdir": "Macrosperma/good_candidates",
        "srnabench_folder_prefix": "Macrosperma_",
        "support_mode": "distinct_libraries",
        **NEMATODE_PROFILE,
    },
    "Sulstoni": {
        "libraries": ["SR0", "SR1", "SR2", "SR3", "SR4", "SR5", "SR6", "SR7"],
        "srnabench_out_subdir": "sRNAtoolboxDB/out/Sulstoni",
        "mirdeep_out_subdir": "Sulstoni/mirdeep_out",
        "scripts_subdir": "Sulstoni/scripts",
        "good_candidates_subdir": "Sulstoni/good_candidates",
        "srnabench_folder_prefix": "Sulstoni_",
        "support_mode": "distinct_libraries",
        **NEMATODE_PROFILE,
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
        **HOFSTENIA_PROFILE,
    },
}


def get_species_config(species, base_path=None, variant=None):
    if species not in SPECIES_CONFIG:
        raise ValueError(
            f"Unknown species '{species}'. "
            f"Choose from: {', '.join(SPECIES_CONFIG.keys())}"
        )
    cfg = SPECIES_CONFIG[species].copy()
    cfg["species"] = species
    cfg["base_path"] = base_path or DEFAULT_BASE_PATH

    if variant == "new_genome":
        if species != "Hofstenia":
            raise ValueError("--variant new_genome is only supported for Hofstenia")
        cfg.update(HOFSTENIA_NEW_GENOME_VARIANT)

    cfg["scripts_dir"] = os.path.join(cfg["base_path"], cfg["scripts_subdir"])
    cfg["good_candidates_dir"] = os.path.join(cfg["base_path"], cfg["good_candidates_subdir"])
    cfg["srnabench_out_dir"] = os.path.join(cfg["base_path"], cfg["srnabench_out_subdir"])
    cfg["mirdeep_out_dir"] = os.path.join(cfg["base_path"], cfg["mirdeep_out_subdir"])
    cfg["seed_path"] = os.path.join(cfg["base_path"], cfg["seed_file"])
    cfg["output_dir"] = os.path.join(
        cfg["base_path"],
        cfg["output_subdir"].format(species=cfg.get("species_label", species)),
    )
    cfg["display_species"] = cfg.get("species_label", species)
    return cfg


def resolve_seed_path(cfg, seed_override=None):
    if seed_override:
        return seed_override
    return cfg["seed_path"]


def load_seed_table(path, encoding="utf-8", sep="\t"):
    import pandas as pd

    if sep == "\t":
        return pd.read_csv(path, sep=sep, encoding=encoding)
    return pd.read_csv(path, encoding=encoding)


def species_uses_blast(cfg):
    return cfg.get("use_blast", False)


def species_uses_mirbase(cfg):
    return cfg.get("use_mirbase_intersects", False)


def species_ziv_unfiltered_only(cfg):
    return cfg.get("ziv_profile") == "unfiltered_only"


def srnabench_folder_name(cfg, library):
    return f"{cfg['srnabench_folder_prefix']}{library}"


def mirdeep_folder_name(cfg, library):
    return library


def dev_condition(cfg, library_name, tool_name):
    """Return grouping key for good_candidates support filtering."""
    lib = library_name.split("_")[-1]
    if cfg["support_mode"] == "dev_condition_replicates":
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
