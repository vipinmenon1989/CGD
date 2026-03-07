#!/usr/bin/env python3
"""
CGD — Comprehensive Guide RNA Design
======================================
Copyright : Vipin Menon & BIG Lab, Hanyang University (HYU)
Author    : Vipin Menon  <a.vipin.menon@gmail.com>
Date      : 18 September 2019  (modernized 2026)

Description
-----------
CGD is a Python-based on-target scoring method for:
  - CRISPRi    (CGDi)
  - CRISPRa    (CGDa)
  - CRISPR-Cas9 canonical   (CGD9)
  - CRISPR-Cas9 non-canonical PAM  (CGD9NG)
  - CRISPR-Cas12a  (CGD12a)

Scores are derived from ENLOR (Elastic Net Logistic Regression) weights
combining free energy, sequence entropy, GC content, and position-dependent
nucleotide features.

Input
-----
A FASTA file with one or more sequences (100 – 10 000 nt each).

Output
------
A tab-separated text file with columns:
    ID  Start  End  Strand  Sequence  <score_column(s)>

Usage
-----
  python CGD.py -a input.fa          # Comprehensive (all systems)
  python CGD.py -b input.fa          # CRISPRi only
  python CGD.py -c input.fa          # CRISPRa only
  python CGD.py -d input.fa          # Cas9 canonical only
  python CGD.py -e input.fa          # Cas12a only
  python CGD.py -f input.fa          # Cas9 non-canonical only
"""

import argparse
import math
import sys
from collections import OrderedDict
from itertools import chain
from typing import Dict, Iterator, List, Optional, Tuple

import numpy as np

try:
    import RNA
except ImportError as exc:
    raise ImportError(
        "ViennaRNA Python bindings are required.\n"
        "Install with:  conda install -c bioconda viennarna\n"
        "           or: pip install ViennaRNA"
    ) from exc

from get_sequence import reverse_complement


# ---------------------------------------------------------------------------
# Type aliases
# ---------------------------------------------------------------------------
# Each guide record: [start, end, strand, sequence, score1, ..., seq_id]
GuideRecord = List


# ===========================================================================
# FASTA parser
# ===========================================================================

def parse_fasta(filepath: str) -> Iterator[Tuple[str, str]]:
    """
    Parse a FASTA file and yield (sequence_id, sequence) pairs.

    Multi-line sequences are concatenated. Each sequence is uppercased and
    stripped of whitespace.

    Args:
        filepath: Path to a FASTA-formatted file.

    Yields:
        Tuples of (sequence_id, full_sequence).

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the file contains no valid FASTA records.
    """
    import os
    if not os.path.isfile(filepath):
        raise FileNotFoundError(f"Input file not found: '{filepath}'")

    current_id: Optional[str] = None
    sequence_parts: List[str] = []
    found_any = False

    with open(filepath, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # Flush previous record
                if current_id is not None:
                    yield current_id, "".join(sequence_parts).upper()
                    found_any = True
                # Start new record — take only the first word as the ID
                current_id = line[1:].split()[0]
                sequence_parts = []
            else:
                sequence_parts.append(line)

    # Flush last record
    if current_id is not None:
        yield current_id, "".join(sequence_parts).upper()
        found_any = True

    if not found_any:
        raise ValueError(f"No valid FASTA records found in '{filepath}'")


# ===========================================================================
# Nucleotide entropy helper
# ===========================================================================

def _sequence_entropy(region: str) -> float:
    """
    Calculate the Shannon entropy (bits) of the nucleotide composition
    of *region*, rounded to one decimal place.

    Args:
        region: A nucleotide string to analyse.

    Returns:
        Rounded Shannon entropy value.
    """
    length = len(region)
    entropy_sum = 0.0
    for nucleotide in ("A", "T", "G", "C"):
        freq = region.count(nucleotide) / float(length)
        if freq > 0:
            entropy_sum -= freq * np.log2(freq)
    return round(entropy_sum, 1)


# ===========================================================================
# Scoring functions
# ===========================================================================

def score_cas12a(sequence: str) -> float:
    """
    Calculate the CINDEL on-target activity score for a Cas12a (Cpf1) guide.

    The model uses free energy, entropy, GC content, dinucleotide composition,
    and 208 position-dependent features trained by ENLOR.

    Args:
        sequence: A 34-bp guide RNA sequence (4 bp upstream + 30 bp target).

    Returns:
        Logistic-transformed activity score in (0, 1).
    """
    # --- Model parameters (position, nucleotide/dinucleotide, weight) ---
    params: List[Tuple[str, int, float]] = [
        ('A',13,0.359277769),('A',28,0.128732671),('A',2,0.089207956),
        ('A',31,0.050508798),('A',7,0.063364247),('A',8,0.193517666),
        ('AA',0,0.074052742),('AA',9,-0.218926976),('AA',10,-0.667432796),
        ('AA',11,-0.521040934),('AA',12,-0.264806378),('AA',19,-0.295058409),
        ('AA',21,-0.17718561),('AA',25,0.186318341),('AA',28,0.176168855),
        ('AA',2,0.11203447),('AC',0,0.119809836),('AC',9,-0.427869295),
        ('AC',13,0.176192535),('AC',18,0.101264344),('AC',1,-0.057025748),
        ('AC',19,0.123838693),('AC',20,0.12864349),('AC',27,-0.115956893),
        ('AC',30,-0.211710227),('AC',31,0.055096611),('AC',32,-0.066466648),
        ('AC',7,0.284818922),('AC',8,0.168980629),('AG',12,0.107615207),
        ('AG',13,0.116406611),('AG',14,0.085343537),('AG',23,-0.052145867),
        ('AG',25,-0.147681142),('AG',32,-0.251946146),('AG',7,0.269547232),
        ('AT',12,-0.055738967),('AT',14,-0.280193723),('AT',16,-0.072299285),
        ('AT',17,-0.188148058),('AT',21,0.13577806),('AT',24,0.071551024),
        ('AT',28,0.121715007),('AT',30,-0.128900537),('C',15,0.187281105),
        ('C',16,0.050394743),('C',17,0.046877049),('C',28,-0.046360015),
        ('C',2,-0.139751449),('C',29,-0.108997003),('C',3,-0.58046435),
        ('C',7,0.282845387),('CA',0,-0.036713083),('CA',9,-0.129619149),
        ('CA',10,-0.228614104),('CA',13,-0.095178226),('CA',1,0.163224633),
        ('CA',19,0.056594182),('CA',20,0.197202638),('CA',21,0.027273837),
        ('CA',22,-0.139136158),('CA',7,0.154463045),('CA',8,-0.112498711),
        ('CC',13,0.129634035),('CC',14,0.053211288),('CC',15,0.056682693),
        ('CC',16,0.099331396),('CC',23,-0.099912427),('CC',28,-0.109688314),
        ('CC',2,-0.207256935),('CC',29,-0.047640371),('CC',32,-0.05591437),
        ('CC',8,0.65410187),('CG',10,-0.193574201),('CG',12,-0.085785037),
        ('CG',13,-0.231165614),('CG',14,-0.11446628),('CG',15,-0.119233092),
        ('CG',17,-0.188450362),('CG',18,-0.090577903),('CG',1,0.755899442),
        ('CG',25,-0.103394296),('CG',26,-0.235079641),('CG',27,-0.196914463),
        ('CG',28,-0.102209147),('CG',2,0.885302791),('CG',32,0.072689139),
        ('CG',8,-0.108629339),('CT',9,0.348496508),('CT',11,-0.224637913),
        ('CT',17,0.054202198),('CT',18,0.081685704),('CT',26,0.033291728),
        ('CT',27,0.119643634),('CT',3,-0.183633657),('G',12,0.027338979),
        ('G',16,-0.053448765),('G',18,-0.073831779),('G',24,-0.137566753),
        ('G',27,-0.350861759),('G',2,0.036045258),('G',3,0.096654517),
        ('G',8,0.826090966),('GA',9,0.359935685),('GA',10,0.078569053),
        ('GA',12,0.169220542),('GA',13,0.058267212),('GA',16,-0.096677117),
        ('GA',25,0.183733491),('GA',30,0.049311266),('GA',31,-0.066188348),
        ('GA',32,0.07418699),('GC',10,-0.372731969),('GC',11,-0.072920116),
        ('GC',12,-0.18347655),('GC',16,0.07684379),('GC',27,-0.207964823),
        ('GC',28,-0.068003699),('GC',30,-0.169119607),('GC',31,-0.084231075),
        ('GG',10,0.061275537),('GG',11,0.363683574),('GG',12,0.107234473),
        ('GG',13,0.2762344),('GG',15,-0.130283364),('GG',16,-0.210606843),
        ('GG',22,-0.278365734),('GG',23,-0.392061429),('GG',24,-0.375657208),
        ('GG',25,-0.276708759),('GG',26,-0.053345569),('GG',32,-0.039499132),
        ('GG',7,0.056485458),('GT',9,0.191596031),('GT',11,0.166350263),
        ('GT',12,0.087029602),('GT',20,0.086971691),('GT',22,0.055831476),
        ('GT',23,0.04992352),('GT',26,-0.036547343),('GT',27,-0.227481683),
        ('GT',28,-0.086145062),('GT',29,-0.316010288),('GT',7,-0.162305425),
        ('GT',8,0.191117867),('T',14,-0.174791539),('T',15,-0.064452789),
        ('T',17,-0.124207175),('T',20,-0.035009888),('T',3,0.233384052),
        ('T',7,-3.464862095),('T',8,-0.184797232),('TA',10,0.311652992),
        ('TA',14,-0.085557183),('TA',16,0.097333539),('TA',18,-0.206798454),
        ('TA',24,0.235579018),('TA',28,0.366518978),('TA',2,-0.457700388),
        ('TA',6,0.064310276),('TA',7,0.460015141),('TC',12,-0.058018063),
        ('TC',16,-0.058701337),('TC',18,0.212081096),('TC',1,-0.108103079),
        ('TC',19,0.134732673),('TC',23,0.045835459),('TC',26,0.245536487),
        ('TC',2,-0.0720811),('TC',29,0.229540677),('TC',30,0.035413033),
        ('TC',6,0.091569704),('TC',7,-0.651599392),('TG',9,0.312514593),
        ('TG',10,0.321561498),('TG',1,0.312174546),('TG',21,0.041875516),
        ('TG',22,0.088872094),('TG',24,0.039480307),('TG',31,-0.099295609),
        ('TT',0,-0.089958105),('TT',9,-0.746490329),('TT',10,-0.500057876),
        ('TT',11,-0.577720541),('TT',12,-0.806473095),('TT',13,-0.752573367),
        ('TT',14,-0.717943612),('TT',15,-0.648672115),('TT',16,-0.353673401),
        ('TT',17,-0.497720724),('TT',1,-0.089333716),('TT',19,-0.557002415),
        ('TT',20,-0.553738685),('TT',21,-0.581588808),('TT',22,-0.442802126),
        ('TT',24,-0.347367156),('TT',25,-0.320309377),('TT',28,0.135942601),
        ('TT',3,0.027108425),('TT',6,-0.304449945),('TT',7,-0.443442269),
        ('TT',8,-1.188734393),
    ]
    intercept     = -1.139975209
    w_free_energy =  0.187763689
    w_entropy     =  1.125236597
    w_AC          =  0.023788317
    w_AG          =  0.013631088
    w_CC          = -0.162320093
    w_CG          =  0.082469431
    w_CT          =  0.039421446
    w_GC          = -0.022381161
    w_GT          =  0.017215267
    w_TA          =  0.084467023
    gc_weight_low  = -0.00000000000499
    gc_weight_high =  0.666115596

    score = intercept

    # Free energy of region [8:31]
    energy_region = sequence[8:31]
    free_energy = round(RNA.fold(energy_region)[-1], 0)
    score += free_energy * w_free_energy

    # Shannon entropy of region [4:27]
    score += _sequence_entropy(sequence[4:27]) * w_entropy

    # Global dinucleotide composition
    score += w_AG * sequence.count("AG")
    score += w_AC * sequence.count("AC")
    score += w_CG * sequence.count("CG")
    score += w_CC * sequence.count("CC")
    score += w_CT * sequence.count("CT")
    score += w_GC * sequence.count("GC")
    score += w_GT * sequence.count("GT")
    score += w_TA * sequence.count("TA")

    # GC-content penalty on the guide region [8:31]
    guide_region = sequence[8:31]
    gc_count = guide_region.count("G") + guide_region.count("C")
    gc_weight = gc_weight_low if gc_count <= 9 else gc_weight_high
    score += abs(9 - gc_count) * gc_weight

    # Position-dependent features
    for nucleotide, position, weight in params:
        if sequence[position: position + len(nucleotide)] == nucleotide:
            score += weight

    return 1.0 / (1.0 + math.exp(-score))


def score_cas9(sequence: str) -> float:
    """
    Calculate the LINDEL on-target score for CRISPR-Cas9 canonical PAM (NGG).

    Args:
        sequence: A 30-bp guide RNA sequence.

    Returns:
        Logistic-transformed activity score in (0, 1).
    """
    parameters: List[Tuple[str, int, float]] = [
        ('A',13,0.1352634445),('A',18,0.2266601183),('A',27,-0.102297985),
        ('T',24,-0.1930543379),('G',4,0.407190242),('G',7,0.2783643584),
        ('G',10,0.0697306019),('G',17,-0.3592921504),('G',18,-0.1287897125),
        ('G',20,0.3691969814),('G',22,0.4827249725),('G',23,1.3789539029),
        ('G',24,0.1259575325),('G',27,-0.2625186279),('G',6,-0.1587906028),
        ('C',21,0.4349944451),('C',27,0.38874478),('AA',18,0.2991883708),
        ('AA',21,-0.4959624082),('TA',17,0.4111542835),('TA',19,-0.2819993431),
        ('TA',20,-1.2192705551),('TA',22,-0.3685329773),('GA',18,0.2622360001),
        ('CA',21,0.3077715591),('CA',22,0.4931654331),('CA',27,0.3573515163),
        ('AT',21,-0.3234310359),('TT',2,0.3494023736),('TT',4,-0.3881133926),
        ('TT',8,-0.2643292347),('TT',13,-0.7733015532),('TT',14,-0.5396509519),
        ('TT',15,-0.2632898469),('TT',19,-0.48701914),('TT',20,-1.6781779772),
        ('TT',21,-1.113814519),('TT',22,-0.9512292771),('TT',23,-0.2276175337),
        ('GT',6,0.3553403317),('GT',17,-0.2247142002),('GT',18,-0.2918636039),
        ('GT',22,1.0425529535),('CT',16,0.2859762363),('CT',17,0.1824673344),
        ('CT',22,-0.7220324868),('CT',23,-0.1886224787),('AG',19,0.301297765),
        ('AG',20,-0.5703135186),('TG',0,-0.1225163563),('TG',9,0.2277042429),
        ('TG',12,-0.2042563995),('TG',16,-0.1867730992),('TG',21,0.3007303605),
        ('GG',12,0.2508979853),('GG',19,-0.3387841934),('CG',16,-0.2579818758),
        ('CG',17,-0.1991746456),('CG',20,0.2107592806),('CG',27,0.2430978627),
        ('AC',14,0.1378218027),('AC',20,0.9705751273),('AC',21,0.7595488933),
        ('TC',5,-0.1842673099),('TC',11,-0.1419670918),('TC',22,-0.6059995886),
        ('GC',9,0.2384841836),('GC',10,0.318694388),('GC',17,-0.2323498858),
        ('GC',18,-0.158530008),('GC',19,0.2060297787),('CC',4,-0.2574206023),
        ('CC',9,-0.1012745036),('CC',23,0.296799014),('CC',28,-0.3137204638),
    ]
    intercept      =  0.675532755
    w_free_energy  =  0.168073832
    w_entropy      =  0.4037850136
    gc_high        =  0.0346278599
    gc_low         = -2.13e-6
    w_T            = -0.0478914162
    w_AC           =  0.0706769127
    w_CG           = -0.0786235386
    w_GT           =  0.0751977394
    w_TT           = -0.2414317576

    score = intercept
    region = sequence[4:24]

    score += round(RNA.fold(region)[-1], 0) * w_free_energy
    score += _sequence_entropy(region) * w_entropy
    score += w_AC * sequence.count("AC")
    score += w_T  * sequence.count("T")
    score += w_CG * sequence.count("CG")
    score += w_GT * sequence.count("GT")
    score += w_TT * sequence.count("TT")

    guide_region = sequence[4:24]
    gc_count = guide_region.count("G") + guide_region.count("C")
    gc_weight = gc_low if gc_count <= 10 else gc_high
    score += abs(10 - gc_count) * gc_weight

    for nucleotide, position, weight in parameters:
        if sequence[position: position + len(nucleotide)] == nucleotide:
            score += weight

    return 1.0 / (1.0 + math.exp(-score))


def score_crispra(sequence: str) -> float:
    """
    Calculate the LINDEL on-target score for CRISPRa (activation).

    Args:
        sequence: A 30-bp guide RNA sequence.

    Returns:
        Logistic-transformed activity score in (0, 1).
    """
    parameters: List[Tuple[str, int, float]] = [
        ("A",20,0.137170389),("T",1,-0.074585823),("T",7,0.111059881),
        ("T",29,0.098637193),("C",7,-0.049469777),("C",9,0.120651133),
        ("C",13,-0.076869685),("C",15,0.095490951),("AA",13,0.441176916),
        ("AA",23,0.151540679),("TA",1,-0.087314896),("TA",4,0.417718773),
        ("TA",7,0.105648267),("TA",8,0.570685229),("TA",10,-0.665609037),
        ("TA",14,-0.185567376),("TA",22,0.869906711),("TA",27,0.108646417),
        ("GA",4,-0.130676946),("GA",8,-0.471092337),("GA",10,-0.190239216),
        ("GA",14,0.158674623),("GA",19,0.251114014),("GA",21,0.286665059),
        ("GA",25,0.414530408),("CA",18,-0.422451619),("CA",28,0.476939986),
        ("AT",8,-0.10205527),("AT",11,-0.135439919),("AT",12,-0.340822076),
        ("AT",15,-0.220870944),("AT",21,0.76460505),("AT",23,-1.125123608),
        ("TT",5,-0.258567901),("TT",8,0.684879543),("TT",11,0.263153909),
        ("TT",13,-0.226042759),("GT",0,-0.145425996),("GT",2,-0.056128866),
        ("GT",6,0.293513389),("GT",7,0.116560365),("GT",17,-0.557781356),
        ("CT",4,-0.143478539),("CT",21,-0.127121675),("CT",22,0.146064932),
        ("CT",23,-0.054563533),("CT",25,0.253908979),("CT",27,-0.160033658),
        ("AG",1,0.181299386),("AG",3,0.088189755),("AG",16,-0.169852859),
        ("TG",8,-0.629010616),("TG",9,-0.532089435),("TG",13,0.116569465),
        ("TG",21,-0.088926423),("GG",0,0.097573308),("GG",4,0.089778821),
        ("GG",13,0.043998425),("GG",17,0.051130828),("GG",20,0.183442263),
        ("CG",12,0.584212426),("AC",5,-0.487374162),("AC",6,-0.191273931),
        ("AC",7,0.849177822),("AC",12,0.277918839),("AC",19,-0.091508321),
        ("TC",4,-0.802078354),("TC",28,-0.144629743),("GC",7,-0.052976793),
        ("GC",15,0.360078914),("GC",22,0.079236972),("CC",0,0.075396858),
        ("CC",5,-0.050765507),("CC",19,0.115084718),("CC",22,-0.089198303),
    ]
    intercept      = -1.671946104
    w_entropy      =  0.4037850136
    gc_high        =  0.240138978
    gc_low         = -0.016037794
    w_free_energy  = -0.000342883
    w_GG           =  0.090650765
    w_AT           = -0.040546895
    w_melt_temp    =  0.005148235

    score = intercept
    region = sequence[4:24]

    score += round(RNA.fold(region)[-1], 0) * w_free_energy
    score += w_AT * sequence.count("AT")
    score += w_GG * sequence.count("GG")

    guide_region = sequence[4:24]
    gc_count = guide_region.count("G") + guide_region.count("C")
    gc_weight = gc_low if gc_count <= 10 else gc_high
    score += abs(10 - gc_count) * gc_weight

    melt_temp = 64.9 + 41.0 * ((guide_region.count("G") + guide_region.count("C") - 16.4) / 20.0)
    score += melt_temp * w_melt_temp

    for nucleotide, position, weight in parameters:
        if sequence[position: position + len(nucleotide)] == nucleotide:
            score += weight

    return 1.0 / (1.0 + math.exp(-score))


def score_crispri(sequence: str) -> float:
    """
    Calculate the LINDEL on-target score for CRISPRi (interference).

    Args:
        sequence: A 30-bp guide RNA sequence.

    Returns:
        Logistic-transformed activity score in (0, 1).
    """
    parameters: List[Tuple[str, int, float]] = [
        ('A',18,0.0815450954),('T',3,0.1251411781),('T',14,-0.05601527),
        ('T',16,-0.1374303579),('T',26,0.6129942086),('G',3,-0.1584387612),
        ('C',22,-0.0687106001),('C',23,-0.1450156038),('C',27,0.1564440206),
        ('C',28,-0.1007689224),('AA',9,-0.1547718709),('AA',11,-0.3012529992),
        ('AA',22,0.1655011092),('TA',8,-0.1481922885),('GA',13,0.1244079638),
        ('GA',18,0.1880705815),('CA',20,-0.3439729695),('CA',21,0.1416105185),
        ('CA',22,0.1522284385),('CA',27,0.1847169166),('CA',28,-0.1975980794),
        ('AT',2,0.3944690605),('AT',8,-0.1972874534),('AT',11,-0.2336100532),
        ('AT',22,-0.2545299014),('AT',27,-0.3396280329),('TT',5,-0.3974796405),
        ('TT',6,-0.0900360596),('TT',7,-0.2787390725),('TT',9,-0.1649851365),
        ('TT',10,-0.2572279591),('TT',11,-0.302645449),('TT',13,-0.3616339187),
        ('TT',14,-0.265153483),('TT',20,-0.293255811),('TT',21,-0.4450813624),
        ('TT',23,-0.2013282932),('GT',0,0.1335707585),('GT',1,-0.1045879319),
        ('GT',7,0.1046889318),('CT',15,-0.1706747529),('CT',23,-0.2440952346),
        ('AG',13,0.0674383094),('AG',21,-0.0548264226),('TG',16,-0.0621481838),
        ('GG',18,-0.3465682018),('GG',19,-0.1770629408),('GG',26,-0.207419321),
        ('GG',27,-0.2534565649),('CG',1,0.1439811675),('CG',7,0.1884952631),
        ('CG',10,0.1179399121),('CG',16,-0.2754010247),('CG',17,-0.1893117011),
        ('AC',15,0.2027064141),('AC',18,0.3311347746),('AC',20,0.2303312739),
        ('AC',28,0.1772425831),('TC',9,-0.1755042319),('GC',1,0.1250831227),
        ('GC',5,0.1045485972),('GC',11,0.0851007713),('GC',19,-0.1945267787),
        ('GC',20,-0.1458484707),('GC',21,-0.5363020499),('GC',22,-0.5256951798),
        ('GC',23,0.1888301689),('GC',25,-0.7042919152),('CC',18,0.1522234528),
        ('CC',22,-0.1439032381),('CC',28,-0.1504334951),
    ]
    intercept      = -1.3484915738
    w_free_energy  =  0.0584654915
    w_entropy      =  0.4056274813
    gc_high        =  0.7542669585
    gc_low         = -0.0065689225
    w_TT           = -0.1045974512
    w_AT           = -0.0957803804
    w_AG           =  0.1051405001
    w_GG           =  0.0459548209
    w_GT           =  0.0463509282
    w_AA           = -0.0437729377
    w_TA           =  0.1324070584

    score = intercept
    region = sequence[4:24]

    score += round(RNA.fold(region)[-1], 0) * w_free_energy
    score += _sequence_entropy(region) * w_entropy
    score += w_AG * sequence.count("AG")
    score += w_AT * sequence.count("AT")
    score += w_GG * sequence.count("GG")
    score += w_TT * sequence.count("TT")
    score += w_TA * sequence.count("TA")
    score += w_AA * sequence.count("AA")
    score += w_GT * sequence.count("GT")

    guide_region = sequence[4:24]
    gc_count = guide_region.count("G") + guide_region.count("C")
    gc_weight = gc_low if gc_count <= 10 else gc_high
    score += abs(10 - gc_count) * gc_weight

    for nucleotide, position, weight in parameters:
        if sequence[position: position + len(nucleotide)] == nucleotide:
            score += weight

    return 1.0 / (1.0 + math.exp(-score))


def score_cas9ng(sequence: str) -> float:
    """
    Calculate the LINDEL on-target score for CRISPR-Cas9 non-canonical PAM (NGA/NGC/NGT).

    Args:
        sequence: A 30-bp guide RNA sequence.

    Returns:
        Logistic-transformed activity score in (0, 1).
    """
    parameters: List[Tuple[str, int, float]] = [
        ('AA',7,-0.166151119),('AA',15,-0.116635449),('AA',16,-0.141978972),
        ('AA',20,-0.46209483),('AA',21,-1.159033893),('AA',23,-0.14092585),
        ('AA',28,-0.140710973),('TA',7,0.242102512),('TA',10,-0.051811653),
        ('TA',16,0.462935638),('TA',17,0.587570687),('TA',21,-0.321626833),
        ('TA',22,0.221199793),('TA',26,-0.308936349),('TA',27,0.228218948),
        ('GA',0,-0.213575878),('GA',6,0.155415387),('GA',11,-0.066555092),
        ('GA',12,0.155827764),('GA',13,0.191500133),('GA',14,0.105536077),
        ('GA',15,0.104278628),('GA',17,-0.134284261),('GA',18,0.144061168),
        ('GA',25,0.111163684),('CA',2,0.038959214),('CA',7,0.126214472),
        ('CA',8,-0.150300426),('CA',17,0.297525515),('CA',19,0.25467499),
        ('CA',20,-0.293330178),('CA',22,0.238690067),('CA',23,-0.21982359),
        ('CA',27,2.054117089),('CA',28,0.335055713),('AT',1,-0.097386892),
        ('AT',2,0.123999347),('AT',8,0.081022603),('AT',10,-0.108078392),
        ('AT',12,-0.079018399),('AT',18,-0.174383894),('AT',19,0.320177513),
        ('TT',6,-0.221409382),('TT',7,-0.157901062),('TT',8,-0.061772525),
        ('TT',10,-0.150198902),('TT',11,-0.448551556),('TT',12,-0.622624981),
        ('TT',13,-0.139416659),('TT',14,-0.077523842),('TT',15,-0.948417398),
        ('TT',16,-0.080861838),('TT',17,-1.046258319),('TT',18,-0.375026297),
        ('TT',20,-1.119104273),('TT',21,-0.166128832),('TT',22,-0.930378086),
        ('GT',21,0.510522765),('GT',22,0.332401486),('GT',23,0.219019899),
        ('GT',25,-0.154193855),('CT',10,0.202473712),('CT',14,-0.158357686),
        ('CT',15,-0.066814565),('CT',16,0.160473088),('CT',22,-0.552932155),
        ('AG',13,0.160018261),('AG',19,-0.108198604),('AG',20,0.144071434),
        ('AG',22,0.25226124),('AG',23,1.063887668),('AG',24,-0.009612897),
        ('AG',28,0.122018273),('TG',5,0.241363252),('TG',9,0.079465841),
        ('TG',10,0.149806008),('TG',12,0.049946857),('TG',13,0.155486915),
        ('TG',14,0.031926277),('TG',18,-0.277324377),('TG',20,0.691550937),
        ('TG',26,-0.08254848),('GG',1,-0.157621236),('GG',6,0.278519061),
        ('GG',9,0.259115271),('GG',23,2.055614334),('CG',2,-0.436407352),
        ('CG',7,-0.046536109),('CG',16,-0.15453076),('CG',17,-0.205884833),
        ('CG',18,0.191552793),('CG',20,-0.245249633),('CG',21,-0.542330163),
        ('CG',27,-0.185411057),('AC',0,-0.085648062),('AC',6,0.148699235),
        ('AC',18,0.105302637),('AC',20,0.175012001),('AC',23,0.264370297),
        ('AC',26,-0.184483722),('TC',0,0.112147888),('TC',11,0.217316946),
        ('TC',13,-0.310569559),('TC',14,-0.121534773),('TC',21,-0.247459342),
        ('TC',22,-1.241334903),('TC',23,-0.628940116),('TC',27,0.536660199),
        ('TC',28,-0.093641005),('GC',1,0.276788432),('GC',4,-0.00869246),
        ('GC',9,0.128560922),('GC',10,0.080325072),('GC',12,-0.219567722),
        ('GC',18,-0.073212082),('GC',19,-0.210102633),('GC',20,-0.457465406),
        ('GC',21,-0.17993102),('GC',22,-1.458409354),('CC',6,-0.110326586),
        ('CC',9,-0.259044689),('CC',16,0.38220164),('CC',17,0.13474772),
        ('CC',18,0.352748051),('CC',19,0.771334455),('CC',21,0.234536174),
        ('CC',26,0.757158333),('CC',28,-0.136943015),
        ('A',10,-0.04673173),('A',20,0.05239417),('A',24,-0.22018074),
        ('A',26,3.09091971),('T',19,-0.36533529),('T',21,-0.35869458),
        ('T',23,-0.35126663),('T',26,-0.8356666),('G',6,0.05688286),
        ('G',17,-0.10189948),('G',20,-0.05641955),('G',23,0.22952347),
        ('G',24,0.90850369),('G',29,0.07510093),('C',5,-0.11228712),
        ('C',17,0.10768524),('C',18,0.05610607),('C',19,0.05474233),
        ('C',20,-0.02714258),('C',21,0.13026953),('C',23,-0.27420231),
    ]
    intercept      = -4.472148961
    w_entropy      =  1.459539385
    w_free_energy  =  0.17806217
    gc_high        =  0.012356139
    gc_low         = -1.05e-13
    w_AG           =  0.09453565
    w_TG           =  0.034280035
    w_TT           = -0.179103124

    score = intercept
    region = sequence[4:24]

    score += round(RNA.fold(region)[-1], 0) * w_free_energy
    score += _sequence_entropy(region) * w_entropy
    score += w_AG * sequence.count("AG")
    score += w_TT * sequence.count("TT")
    score += w_TG * sequence.count("TA")   # Note: original uses TG weight on TA count

    guide_region = sequence[4:24]
    gc_count = guide_region.count("G") + guide_region.count("C")
    gc_weight = gc_low if gc_count <= 10 else gc_high
    score += abs(10 - gc_count) * gc_weight

    for nucleotide, position, weight in parameters:
        if sequence[position: position + len(nucleotide)] == nucleotide:
            score += weight

    return 1.0 / (1.0 + math.exp(-score))


# ===========================================================================
# PAM scanning helpers
# ===========================================================================

def _scan_forward(
    sequence: str,
    pam_motif: str,
    upstream: int,
    downstream: int,
) -> Iterator[Tuple[int, int, str]]:
    """
    Yield (start, end, window) for every occurrence of *pam_motif* on the
    forward strand where a full upstream+downstream window fits.

    Args:
        sequence:   Full genomic sequence (uppercase).
        pam_motif:  The PAM dinucleotide/trinucleotide to search for.
        upstream:   Bases to include before the PAM position.
        downstream: Bases to include after the PAM position start.

    Yields:
        Tuples of (window_start, window_end, window_sequence).
    """
    search_pos = 0
    while True:
        pos = sequence.find(pam_motif, search_pos)
        if pos == -1:
            break
        if upstream < pos < len(sequence) - downstream:
            start = pos - upstream
            end   = pos + downstream
            yield start, end, sequence[start:end]
        search_pos = pos + 1


def _scan_reverse(
    sequence: str,
    rc_pam_motif: str,
    upstream: int,
    downstream: int,
) -> Iterator[Tuple[int, int, str]]:
    """
    Yield (start, end, rc_window) for every occurrence of the reverse-complement
    PAM motif (*rc_pam_motif*) on the forward strand.

    Args:
        sequence:      Full genomic sequence (uppercase).
        rc_pam_motif:  The reverse-complement of the PAM motif.
        upstream:      Bases upstream of the RC-PAM position start.
        downstream:    Bases downstream of the RC-PAM position start.

    Yields:
        Tuples of (window_start, window_end, reverse_complement_window).
    """
    search_pos = 0
    while True:
        pos = sequence.find(rc_pam_motif, search_pos)
        if pos == -1:
            break
        if upstream < pos < len(sequence) - downstream:
            start = pos - upstream
            end   = pos + downstream
            window = sequence[start:end]
            yield start, end, reverse_complement(window)
        search_pos = pos + 1


# ===========================================================================
# Pipeline functions
# ===========================================================================

def _write_tsv(
    filepath: str,
    header: List[str],
    records: Dict[str, GuideRecord],
    sort_col: int,
) -> None:
    """
    Write guide records to a tab-separated file, sorted by *sort_col* descending.

    Args:
        filepath:  Destination file path.
        header:    Column header list.
        records:   Dict of guide_id → record list.
        sort_col:  Index in the record list to sort by.
    """
    sorted_records = OrderedDict(
        sorted(records.items(), key=lambda item: item[1][sort_col], reverse=True)
    )
    with open(filepath, "w", encoding="utf-8") as out:
        out.write("\t".join(header) + "\n")
        for record in sorted_records.values():
            out.write("\t".join(str(v) for v in record) + "\n")


def run_comprehensive(fasta_path: str) -> None:
    """
    Score all guide RNA candidates (CRISPRi, CRISPRa, Cas9, Cas9NG, Cas12a)
    from a FASTA file and write results to ``CGD.txt``.

    Args:
        fasta_path: Path to input FASTA file (sequences must be 100–10 000 nt).
    """
    records: Dict[str, GuideRecord] = {}

    for seq_id, sequence in parse_fasta(fasta_path):
        if not (100 <= len(sequence) <= 10000):
            print(f"WARNING: Sequence '{seq_id}' length {len(sequence)} is out of range (100–10 000 nt). Skipping.")
            continue

        # --- Cas9 canonical (NGG / NCC reverse) ---
        for start, end, window in _scan_forward(sequence, "GG", 25, 5):
            guide_id = f"gRNAc_{start}"
            records[guide_id] = [
                seq_id, start, end, "+", window,
                round(score_crispri(window), 2),
                round(score_crispra(window), 2),
                round(score_cas9(window), 2),
                0.0, 0.0,
            ]
        for start, end, rc_window in _scan_reverse(sequence, "CC", 3, 27):
            guide_id = f"gRNAcn_{start}"
            records[guide_id] = [
                seq_id, start, end, "-", rc_window,
                round(score_crispri(rc_window), 2),
                round(score_crispra(rc_window), 2),
                round(score_cas9(rc_window), 2),
                0.0, 0.0,
            ]

        # --- Cas9 non-canonical NGA ---
        for start, end, window in _scan_forward(sequence, "GA", 25, 5):
            guide_id = f"gRNAnga_{start}"
            records[guide_id] = [seq_id, start, end, "+", window, 0.0, 0.0, 0.0, round(score_cas9ng(window), 2), 0.0]
        for start, end, rc_window in _scan_reverse(sequence, "TC", 3, 27):
            guide_id = f"gRNAngac_{start}"
            records[guide_id] = [seq_id, start, end, "-", rc_window, 0.0, 0.0, 0.0, round(score_cas9ng(rc_window), 2), 0.0]

        # --- Cas9 non-canonical NGC ---
        for start, end, window in _scan_forward(sequence, "GC", 25, 5):
            guide_id = f"gRNAngc_{start}"
            records[guide_id] = [seq_id, start, end, "+", window, 0.0, 0.0, 0.0, round(score_cas9ng(window), 2), 0.0]
        for start, end, rc_window in _scan_reverse(sequence, "GC", 3, 27):
            guide_id = f"gRNAngcr_{start}"
            records[guide_id] = [seq_id, start, end, "-", rc_window, 0.0, 0.0, 0.0, round(score_cas9ng(rc_window), 2), 0.0]

        # --- Cas9 non-canonical NGT ---
        for start, end, window in _scan_forward(sequence, "GT", 25, 5):
            guide_id = f"gRNAngt_{start}"
            records[guide_id] = [seq_id, start, end, "+", window, 0.0, 0.0, 0.0, round(score_cas9ng(window), 2), 0.0]
        for start, end, rc_window in _scan_reverse(sequence, "AC", 3, 27):
            guide_id = f"gRNAngtr_{start}"
            records[guide_id] = [seq_id, start, end, "-", rc_window, 0.0, 0.0, 0.0, round(score_cas9ng(rc_window), 2), 0.0]

        # --- Cas12a (TTT / AAA reverse) ---
        for start, end, window in _scan_forward(sequence, "TTT", 4, 30):
            guide_id = f"gRNAcas_{start}"
            records[guide_id] = [seq_id, start, end, "+", window, 0.0, 0.0, 0.0, 0.0, round(score_cas12a(window), 2)]
        for start, end, rc_window in _scan_reverse(sequence, "AAA", 27, 7):
            guide_id = f"gRNAcasr_{start}"
            records[guide_id] = [seq_id, start, end, "-", rc_window, 0.0, 0.0, 0.0, 0.0, round(score_cas12a(rc_window), 2)]

    header = ["ID", "Start", "End", "Strand", "Sequence", "CGDi", "CGDa", "CGD9", "CGDNG", "CGD12a"]
    _write_tsv("CGD.txt", header, records, sort_col=7)
    print(f"Comprehensive scoring complete → CGD.txt  ({len(records)} guides)")


def run_cgdi(fasta_path: str) -> None:
    """
    Score CRISPRi guide RNAs from a FASTA file. Output written to ``CGDi.txt``.

    Args:
        fasta_path: Path to input FASTA file (sequences must be 100–10 000 nt).
    """
    records: Dict[str, GuideRecord] = {}

    for seq_id, sequence in parse_fasta(fasta_path):
        if not (100 <= len(sequence) <= 10000):
            print(f"WARNING: '{seq_id}' skipped (length out of range).")
            continue
        for start, end, window in _scan_forward(sequence, "GG", 25, 5):
            records[f"gRNAi_{start}"] = [seq_id, start, end, "+", window, round(score_crispri(window), 2)]
        for start, end, rc_window in _scan_reverse(sequence, "CC", 3, 27):
            records[f"gRNAin_{start}"] = [seq_id, start, end, "-", rc_window, round(score_crispri(rc_window), 2)]

    _write_tsv("CGDi.txt", ["ID", "Start", "End", "Strand", "Sequence", "CGDi"], records, sort_col=5)
    print(f"CRISPRi scoring complete → CGDi.txt  ({len(records)} guides)")


def run_cgda(fasta_path: str) -> None:
    """
    Score CRISPRa guide RNAs from a FASTA file. Output written to ``CGDa.txt``.

    Args:
        fasta_path: Path to input FASTA file (sequences must be 100–10 000 nt).
    """
    records: Dict[str, GuideRecord] = {}

    for seq_id, sequence in parse_fasta(fasta_path):
        if not (100 <= len(sequence) <= 10000):
            print(f"WARNING: '{seq_id}' skipped (length out of range).")
            continue
        for start, end, window in _scan_forward(sequence, "GG", 25, 5):
            records[f"gRNAa_{start}"] = [seq_id, start, end, "+", window, round(score_crispra(window), 2)]
        for start, end, rc_window in _scan_reverse(sequence, "CC", 3, 27):
            records[f"gRNAan_{start}"] = [seq_id, start, end, "-", rc_window, round(score_crispra(rc_window), 2)]

    _write_tsv("CGDa.txt", ["ID", "Start", "End", "Strand", "Sequence", "CGDa"], records, sort_col=5)
    print(f"CRISPRa scoring complete → CGDa.txt  ({len(records)} guides)")


def run_cgd9(fasta_path: str) -> None:
    """
    Score CRISPR-Cas9 guide RNAs from a FASTA file. Output written to ``CGD9.txt``.

    Args:
        fasta_path: Path to input FASTA file (sequences must be 100–10 000 nt).
    """
    records: Dict[str, GuideRecord] = {}

    for seq_id, sequence in parse_fasta(fasta_path):
        if not (100 <= len(sequence) <= 10000):
            print(f"WARNING: '{seq_id}' skipped (length out of range).")
            continue
        for start, end, window in _scan_forward(sequence, "GG", 25, 5):
            records[f"gRNA9_{start}"] = [seq_id, start, end, "+", window, round(score_cas9(window), 2)]
        for start, end, rc_window in _scan_reverse(sequence, "CC", 3, 27):
            records[f"gRNA9n_{start}"] = [seq_id, start, end, "-", rc_window, round(score_cas9(rc_window), 2)]

    _write_tsv("CGD9.txt", ["ID", "Start", "End", "Strand", "Sequence", "CGD9"], records, sort_col=5)
    print(f"Cas9 scoring complete → CGD9.txt  ({len(records)} guides)")


def run_cgd12a(fasta_path: str) -> None:
    """
    Score CRISPR-Cas12a guide RNAs from a FASTA file. Output written to ``CGD12a.txt``.

    Args:
        fasta_path: Path to input FASTA file (sequences must be 100–10 000 nt).
    """
    records: Dict[str, GuideRecord] = {}

    for seq_id, sequence in parse_fasta(fasta_path):
        if not (100 <= len(sequence) <= 10000):
            print(f"WARNING: '{seq_id}' skipped (length out of range).")
            continue
        for start, end, window in _scan_forward(sequence, "TTT", 4, 30):
            records[f"gRNAcas_{start}"] = [seq_id, start, end, "+", window, round(score_cas12a(window), 2)]
        for start, end, rc_window in _scan_reverse(sequence, "AAA", 27, 7):
            records[f"gRNAcasr_{start}"] = [seq_id, start, end, "-", rc_window, round(score_cas12a(rc_window), 2)]

    _write_tsv("CGD12a.txt", ["ID", "Start", "End", "Strand", "Sequence", "CGD12a"], records, sort_col=5)
    print(f"Cas12a scoring complete → CGD12a.txt  ({len(records)} guides)")


def run_cgd9ng(fasta_path: str) -> None:
    """
    Score CRISPR-Cas9 non-canonical PAM guide RNAs (NGA, NGC, NGT) from a FASTA
    file. Output written to ``CGD9NG.txt``.

    Args:
        fasta_path: Path to input FASTA file (sequences must be 100–10 000 nt).
    """
    records: Dict[str, GuideRecord] = {}

    for seq_id, sequence in parse_fasta(fasta_path):
        if not (100 <= len(sequence) <= 10000):
            print(f"WARNING: '{seq_id}' skipped (length out of range).")
            continue
        for pam, rc_pam in [("GA", "TC"), ("GC", "GC"), ("GT", "AC")]:
            for start, end, window in _scan_forward(sequence, pam, 25, 5):
                records[f"gRNAng{pam}_{start}"] = [seq_id, start, end, "+", window, round(score_cas9ng(window), 2)]
            for start, end, rc_window in _scan_reverse(sequence, rc_pam, 3, 27):
                records[f"gRNAngr{pam}_{start}"] = [seq_id, start, end, "-", rc_window, round(score_cas9ng(rc_window), 2)]

    _write_tsv("CGD9NG.txt", ["ID", "Start", "End", "Strand", "Sequence", "CGD9NG"], records, sort_col=5)
    print(f"Cas9-NG scoring complete → CGD9NG.txt  ({len(records)} guides)")


# ===========================================================================
# CLI
# ===========================================================================

def _build_parser() -> argparse.ArgumentParser:
    """Build and return the CGD argument parser."""
    parser = argparse.ArgumentParser(
        prog="CGD",
        description=(
            "CGD — Comprehensive Guide RNA Design\n"
            "On-target scoring for CRISPRi, CRISPRa, Cas9, Cas9-NG, and Cas12a."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python CGD.py -a input.fa          # All systems (comprehensive)\n"
            "  python CGD.py -b input.fa          # CRISPRi only\n"
            "  python CGD.py -c input.fa          # CRISPRa only\n"
            "  python CGD.py -d input.fa          # Cas9 canonical only\n"
            "  python CGD.py -e input.fa          # Cas12a only\n"
            "  python CGD.py -f input.fa          # Cas9 non-canonical only\n"
        ),
    )
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("-a", metavar="FASTA", help="Comprehensive score (all CRISPR systems)")
    mode.add_argument("-b", metavar="FASTA", help="CRISPRi score (CGDi)")
    mode.add_argument("-c", metavar="FASTA", help="CRISPRa score (CGDa)")
    mode.add_argument("-d", metavar="FASTA", help="CRISPR-Cas9 canonical score (CGD9)")
    mode.add_argument("-e", metavar="FASTA", help="CRISPR-Cas12a score (CGD12a)")
    mode.add_argument("-f", metavar="FASTA", help="CRISPR-Cas9 non-canonical score (CGD9NG)")
    return parser


def main() -> None:
    """Entry point for the CGD command-line tool."""
    parser = _build_parser()
    args = parser.parse_args()

    try:
        if args.a:
            run_comprehensive(args.a)
        elif args.b:
            run_cgdi(args.b)
        elif args.c:
            run_cgda(args.c)
        elif args.d:
            run_cgd9(args.d)
        elif args.e:
            run_cgd12a(args.e)
        elif args.f:
            run_cgd9ng(args.f)
    except (FileNotFoundError, ValueError) as err:
        print(f"ERROR: {err}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
