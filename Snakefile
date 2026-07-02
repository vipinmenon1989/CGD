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

configfile: "config/config.yaml"

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
    params:
        flag=lambda wc: MODE_FLAG[SUFFIX_TO_MODE[wc.suffix]],
    conda:
        "envs/environment.yaml"
    shell:
        "python {workflow.basedir}/CGD.py -{params.flag} {input.fasta} "
        "-o {output} > {log} 2>&1"
