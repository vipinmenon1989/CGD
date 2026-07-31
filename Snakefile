"""
CGD Snakemake workflow
=======================
Runs the CGD on-target guide-RNA scoring tool (CGD.py) over one or more
FASTA inputs, for one or more CRISPR-system scoring modes, and writes
tab-separated results to `results/<sample>/<SCORE>.txt`.

Configure inputs/modes in config/config.yaml, then run:

    snakemake --cores 1 --use-conda

or, if dependencies (numpy, ViennaRNA) are already installed in the
active environment:

    snakemake --cores 1
"""

from pathlib import Path

from snakemake.utils import validate


configfile: "config/config.yaml"
validate(config, workflow.basedir + "/config/config.schema.yaml")

# Mode name -> CGD.py CLI flag
MODE_FLAG = {
    "comprehensive": "a",
    "crispri": "b",
    "crispra": "c",
    "cas9": "d",
    "cas12a": "e",
    "cas9_ng": "f",
}

# Mode name -> output score-column prefix (matches CGD.py's own file naming)
MODE_SUFFIX = {
    "comprehensive": "CGD",
    "crispri": "CGDi",
    "crispra": "CGDa",
    "cas9": "CGD9",
    "cas12a": "CGD12a",
    "cas9_ng": "CGD9NG",
}
SUFFIX_TO_MODE = {v: k for k, v in MODE_SUFFIX.items()}

SAMPLES = config["samples"]
MODES = config.get("modes", list(MODE_FLAG))
OUTDIR = config.get("outdir", "results")
CGD_SCRIPT = str(Path(workflow.basedir) / "CGD.py")

for m in MODES:
    if m not in MODE_FLAG:
        raise ValueError(
            f"Unknown mode '{m}' in config.yaml. "
            f"Valid modes: {sorted(MODE_FLAG)}"
        )


rule all:
    input:
        expand(
            OUTDIR + "/{sample}/{suffix}.txt",
            sample=SAMPLES.keys(),
            suffix=[MODE_SUFFIX[m] for m in MODES],
        ),


rule cgd_score:
    """Score every guide RNA candidate for one sample/mode combination."""
    input:
        fasta=lambda wc: SAMPLES[wc.sample],
    output:
        OUTDIR + "/{sample}/{suffix}.txt",
    log:
        OUTDIR + "/logs/{sample}.{suffix}.log",
    threads: 1
    resources:
        mem_mb=2000,
        runtime=30,
    params:
        flag=lambda wc: MODE_FLAG[SUFFIX_TO_MODE[wc.suffix]],
        script=CGD_SCRIPT,
    conda:
        "envs/environment.yaml"
    shell:
        "python {params.script:q} -{params.flag} {input.fasta:q} "
        "-o {output:q} > {log:q} 2>&1"
