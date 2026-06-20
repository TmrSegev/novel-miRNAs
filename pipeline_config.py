"""Shared configuration for per-library miRNA discovery pipelines."""

import os

DEFAULT_BASE_PATH = "/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq"
DEFAULT_NCRNA_DIR = "/mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/ncRNAs_Caenorhabditis"

# Global filter defaults (all species; overridable via CLI).
FILTER_MC_DEFAULT = 10
EXCLUDE_C_DEFAULT = 100
LOW_SCORE_TOTAL_DEFAULT = 100

# miRBase-based nematode miRNA seeds (Elegans, Macrosperma, Sulstoni).
SEED_FILE_NEMATODE = "mirbase_data/Seeds.txt"
# General miRs of all species (mirGeneDB family table; Hofstenia).
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

NEW_GENOME_SUFFIX = "_newGenome"

NEW_GENOME_OVERRIDES = {
    "Hofstenia": {
        "genome_fasta_subpath": "Hofstenia_newGenome/sRNA_PBonly/hofPB_v6.FINAL.fa",
    },
    "Macrosperma": {
        "genome_fasta_subpath": (
            "Macrosperma_newGenome/genome/"
            "CMACR.caenorhabditis_macrosperma_JU2083_v2.scaffolds.fna"
        ),
    },
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


def resolve_species_and_variant(species, variant=None):
    """Resolve Species_newGenome alias to (base_species, new_genome)."""
    if species.endswith(NEW_GENOME_SUFFIX):
        base = species[: -len(NEW_GENOME_SUFFIX)]
        if base in SPECIES_CONFIG:
            if variant and variant != "new_genome":
                raise ValueError(
                    f"Species alias '{species}' implies variant new_genome; got '{variant}'"
                )
            return base, "new_genome"
    if species not in SPECIES_CONFIG:
        valid = ", ".join(sorted(SPECIES_CONFIG.keys()))
        raise ValueError(
            f"Unknown species '{species}'. Choose from: {valid} "
            f"(or {{Species}}{NEW_GENOME_SUFFIX} for alternate assemblies)"
        )
    return species, variant


def build_new_genome_variant(species):
    track = f"{species}{NEW_GENOME_SUFFIX}"
    return {
        "mirdeep_out_subdir": f"{track}/mirdeep_out",
        "scripts_subdir": f"{track}/scripts",
        "good_candidates_subdir": f"{track}/good_candidates",
        "srnabench_out_subdir": f"sRNAtoolboxDB/out/{track}",
        "output_subdir": f"RNAcentral/miRNAs/{track}",
        "species_label": species,
        "variant_track": track,
        "variant_root_subdir": track,
    }


def get_species_config(species, base_path=None, variant=None):
    species, variant = resolve_species_and_variant(species, variant)
    cfg = SPECIES_CONFIG[species].copy()
    cfg["species"] = species
    cfg["base_path"] = base_path or DEFAULT_BASE_PATH
    cfg["variant"] = variant

    if variant == "new_genome":
        cfg.update(build_new_genome_variant(species))
        overrides = NEW_GENOME_OVERRIDES.get(species, {})
        cfg.update(overrides)
        if overrides.get("genome_fasta_subpath"):
            cfg["genome_fasta"] = os.path.join(
                cfg["base_path"], overrides["genome_fasta_subpath"]
            )

    cfg["scripts_dir"] = os.path.join(cfg["base_path"], cfg["scripts_subdir"])
    cfg["good_candidates_dir"] = os.path.join(cfg["base_path"], cfg["good_candidates_subdir"])
    cfg["srnabench_out_dir"] = os.path.join(cfg["base_path"], cfg["srnabench_out_subdir"])
    cfg["mirdeep_out_dir"] = os.path.join(cfg["base_path"], cfg["mirdeep_out_subdir"])
    cfg["seed_path"] = os.path.join(cfg["base_path"], cfg["seed_file"])
    cfg["output_dir"] = os.path.join(
        cfg["base_path"],
        cfg["output_subdir"].format(species=cfg.get("species_label", species)),
    )
    if variant == "new_genome":
        cfg["variant_root_dir"] = os.path.join(cfg["base_path"], cfg["variant_root_subdir"])
        cfg["ziv_output_dir"] = os.path.join(cfg["base_path"], "Ziv_Features")
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
