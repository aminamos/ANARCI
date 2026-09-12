"""End-to-end tests for the anarci() entry point.

These mirror the higher-level tests we run internally but exercise only this
fork's public surface: the ``anarci()`` tuple return, ``run_germline_assignment``
(via ``assign_germline=True``) and the shipped IMGT HMM/germline data. No typed
datamodel wrappers are used here so the tests stay self-contained.

Germline identities are pinned against the corrected IMGT germlines shipped in
this fork, so they double as a regression guard on the germline database.
"""

import os

import pyhmmer
import pytest

from anarci import anarci, validate_sequence
from anarci.anarci import HMM_path, all_germlines, all_species

# A paired antibody (12e8), a single-chain Fv with both a light and a heavy
# domain, and a non-antibody (lysozyme) that must yield no domains.
SEQ_12E8_H = "EVQLQQSGAEVVRSGASVKLSCTASGFNIKDYYIHWVKQRPEKGLEWIGWIDPEIGDTEYVPKFQGKATMTADTSSNTAYLQLSSLTSEDTAVYYCNAGHDYDRGRFPYWGQGTLVTVSAAKTTPPSVYPLAP"
SEQ_12E8_L = "DIVMTQSQKFMSTSVGDRVSITCKASQNVGTAVAWYQQKPGQSPKLMIYSASNRYTGVPDRFTGSGSGTDFTLTISNMQSEDLADYFCQQYSSYPLTFGAGTKLELKRADAAPTVSIFPPSSEQLTSGGASV"
SEQ_SCFV_A = (
    "DIQMTQSPSSLSASVGDRVTITCRTSGNIHNYLTWYQQKPGKAPQLLIYNAKTLADGVPSRFSGSGSGTQFTLTISSLQPEDFANYYCQHFWSLPFTFGQGTKVEIKRTG"
    "GGGSGGGGSGGGGSGGGGSEVQLVESGGGLVQPGGSLRLSCAASGFDFSRYDMSWVRQAPGKRLEWVAYISSGGGSTYFPDTVKGRFTISRDNAKNTLYLQMNSLRAEDT"
    "AVYYCARQNKKLTWFDYWGQGTLVTVSSHHHHHH"
)
SEQ_LYSOZYME_A = "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL"

SEQUENCES = [
    ("12e8:H", SEQ_12E8_H),
    ("12e8:L", SEQ_12E8_L),
    ("scfv:A", SEQ_SCFV_A),
    ("lysozyme:A", SEQ_LYSOZYME_A),
]

# A VHH (nanobody) with a leading signal peptide and trailing constant region,
# and an unusually framed sequence used to exercise subsequence renumbering.
SEQ_8QOT_H = (
    "MKKNIAFLLASMFVFSIATNAYAEISEVQLVESGGGLVQPGGSLRLSCAASGFNFSYYSIHWVRQAPGKGLEWVAYISSSSSYTSYADSVKGRFTISADTSKNTAYLQMN"
    "SLRAEDTAVYYCARGYQYWQYHASWYWNGGLDYWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLS"
    "SVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHT"
)
SEQ_WEIRD_H = (
    "YVSPLSVALGETARISCGRQALGSRAVQWYQHKPGQAPILLIYNNQDRPSGIPERFSGTPDINFGTTATLTISGVEVGDEADYYCHMWDSRSGFSWSFGGATRLTVLSQP"
    "KAAPSVTLFPPSSEELQANKATLVCLISDFYPGAVTVAWKADSSPVKAGVETTTPSKQSNNKYAASSYLSLTPEQWKSHKSYSCQVTHEGSTVEKTVAPT"
)


def _germlines(detail):
    """Flatten a domain's germline assignment into a comparable dict."""
    (v_species, v_gene), v_identity = detail["germlines"]["v_gene"]
    (j_species, j_gene), j_identity = detail["germlines"]["j_gene"]
    return {
        "v_species": v_species,
        "v_gene": v_gene,
        "v_identity": v_identity,
        "j_species": j_species,
        "j_gene": j_gene,
        "j_identity": j_identity,
    }


def _number_to_residue_index(numbered_domain):
    """Map IMGT numbers (e.g. "111A") to their index in the original sequence."""
    numbering, start, _end = numbered_domain
    mapping = {}
    residue_index = start
    for (position, insertion), amino_acid in numbering:
        if amino_acid != "-":
            mapping[f"{position}{insertion}".strip()] = residue_index
            residue_index += 1
    return mapping


def _hmm_name(hmm):
    return hmm.name.decode() if isinstance(hmm.name, bytes) else hmm.name


def test_detect_multiple_ig_domains():
    numbered, alignment_details, hit_tables = anarci(SEQUENCES, scheme="imgt", output=False, assign_germline=True)

    assert len(numbered) == len(alignment_details) == len(hit_tables) == len(SEQUENCES)

    # 12e8:H - one heavy domain.
    assert len(numbered[0]) == 1
    assert alignment_details[0][0]["chain_type"] == "H"
    assert _germlines(alignment_details[0][0]) == {
        "v_species": "mouse",
        "v_gene": "IGHV14-4*02",
        "v_identity": pytest.approx(0.9489795918367347),
        "j_species": "mouse",
        "j_gene": "IGHJ3*01",
        "j_identity": pytest.approx(0.9230769230769231),
    }

    # 12e8:L - one light (kappa) domain.
    assert len(numbered[1]) == 1
    assert alignment_details[1][0]["chain_type"] in {"K", "L"}
    assert _germlines(alignment_details[1][0]) == {
        "v_species": "mouse",
        "v_gene": "IGKV6-13*01",
        "v_identity": pytest.approx(0.9891304347826086),
        "j_species": "mouse",
        "j_gene": "IGKJ5*01",
        "j_identity": pytest.approx(1.0),
    }

    # scfv:A - light domain then heavy domain.
    assert len(numbered[2]) == 2
    assert alignment_details[2][0]["chain_type"] in {"K", "L"}
    assert alignment_details[2][1]["chain_type"] == "H"
    assert _germlines(alignment_details[2][0]) == {
        "v_species": "mouse",
        "v_gene": "IGKV12-41*01",
        "v_identity": pytest.approx(0.8478260869565217),
        "j_species": "mouse",
        "j_gene": "IGKJ4*02",
        "j_identity": pytest.approx(0.8333333333333334),
    }
    assert _germlines(alignment_details[2][1]) == {
        "v_species": "mouse",
        "v_gene": "IGHV5-12-1*01",
        "v_identity": pytest.approx(0.8877551020408163),
        "j_species": "mouse",
        "j_gene": "IGHJ4*01",
        "j_identity": pytest.approx(0.9230769230769231),
    }

    # lysozyme:A - not an antibody, so no domains.
    assert numbered[3] is None
    assert alignment_details[3] is None


def test_domain_sequences_are_contiguous_slices():
    numbered, _alignment_details, _hit_tables = anarci(SEQUENCES, scheme="imgt", output=False)
    for (_name, sequence), domains in zip(SEQUENCES, numbered):
        if domains is None:
            continue
        for numbering, start, end in domains:
            reproduced = "".join(amino_acid for _, amino_acid in numbering if amino_acid != "-")
            assert reproduced == sequence[start : end + 1]


def test_rhesus_imgt_hmms_have_cys_anchor_columns():
    """Corrected rhesus HMMs must place the conserved cysteines at IMGT 23 and 104."""
    with pyhmmer.plan7.HMMFile(os.path.join(HMM_path, "ALL.hmm")) as hmm_file:
        rhesus_hmms = {_hmm_name(hmm): hmm for hmm in hmm_file if _hmm_name(hmm) in {"rhesus_H", "rhesus_K"}}

    assert set(rhesus_hmms) == {"rhesus_H", "rhesus_K"}
    for hmm in rhesus_hmms.values():
        consensus = hmm.consensus
        assert [i + 1 for i, residue in enumerate(consensus) if residue.upper() == "C"] == [23, 104]


def test_rhesus_v_germlines_have_cys_anchor_columns():
    """Every rhesus V germline must keep the conserved cysteines at columns 23 and 104."""
    for _chain_type, species_germlines in all_germlines["V"].items():
        if "rhesus" not in species_germlines:
            continue
        for sequence in species_germlines["rhesus"].values():
            assert sequence[22] == "C"
            assert sequence[103] == "C"


def test_light_chain_conserved_cys_numbering():
    """A human kappa chain must number its conserved cysteines at IMGT 23 and 104.

    Before the germline/HMM reframe this human kappa best-hit the rhesus HMM and
    the anchors drifted; with the corrected data it hits human_K and the
    cysteines land on their canonical columns.
    """
    sequence = "EIVLTQSPGTLSLSPGERATLSCRASQSVSSSFLAWYQQKPGQAPRLLIYYASSRATGIPDRFSGSGSGTDFTLTISRLEPEDFAVYYCQQTGRIPPTFGQGTKVEIK"
    numbered, alignment_details, _hit_tables = anarci([("7ah1_L", sequence)], scheme="imgt", output=False, assign_germline=True)

    detail = alignment_details[0][0]
    assert detail["chain_type"] == "K"
    assert detail["id"] == "human_K"
    assert _germlines(detail) == {
        "v_species": "human",
        "v_gene": "IGKV3-20*01",
        "v_identity": pytest.approx(0.956989247311828),
        "j_species": "human",
        "j_gene": "IGKJ1*01",
        "j_identity": pytest.approx(0.9166666666666666),
    }

    number_to_residue_index = _number_to_residue_index(numbered[0][0])
    assert sequence[number_to_residue_index["23"]] == "C"
    assert sequence[number_to_residue_index["104"]] == "C"


def test_vhh_nanobody():
    """A VHH should number as a single heavy domain with a significant hit."""
    numbered, alignment_details, _hit_tables = anarci([("8qot_H", SEQ_8QOT_H)], scheme="imgt", output=False, assign_germline=True)

    assert len(numbered[0]) == 1
    detail = alignment_details[0][0]
    assert detail["chain_type"] == "H"
    assert detail["evalue"] > 0
    assert detail["bitscore"] > 80


def test_subsequence_numbering_is_consistent():
    """Renumbering an extracted domain must reproduce the parent's residue identities."""
    numbered, _alignment_details, _hit_tables = anarci([("weird_H", SEQ_WEIRD_H)], scheme="imgt", output=False)
    parent_domain = numbered[0][0]
    parent_index = _number_to_residue_index(parent_domain)

    _numbering, start, end = parent_domain
    domain_sequence = SEQ_WEIRD_H[start : end + 1]
    sub_numbered, _sub_details, _sub_hits = anarci([("subsequence", domain_sequence)], scheme="imgt", output=False)
    sub_index = _number_to_residue_index(sub_numbered[0][0])

    for imgt_number, sub_residue_index in sub_index.items():
        if imgt_number not in parent_index:
            continue
        assert SEQ_WEIRD_H[parent_index[imgt_number]] == domain_sequence[sub_residue_index]

    # Conserved FWR1/FWR3 cysteines.
    assert SEQ_WEIRD_H[parent_index["23"]] == "C"
    assert SEQ_WEIRD_H[parent_index["104"]] == "C"


def test_error_handling():
    """validate_sequence rejects non-amino-acid characters and over-long input."""
    assert validate_sequence("EVQLQQSGAEVVRSGASVKLSCTASGFNIK") is True

    for bad_sequence in ("EVQLQQSGAEVVRSGA123VKLSCTASGFNIK", "EVQLQQSGAEVVRSGASVKLSCTASGFNIK@#$", "EVQLQQSGAEVVRSGASVKLSCTASGFNIKX"):
        with pytest.raises(AssertionError):
            validate_sequence(bad_sequence)

    with pytest.raises(AssertionError):
        validate_sequence("A" * 11000)


def test_all_species_matches_germline_database():
    """all_species must mirror the species present in the germline V database.

    If this fails, a species was added to or removed from the shipped germline
    database and downstream consumers relying on all_species need updating.
    """
    germline_species = set()
    for species_germlines in all_germlines["V"].values():
        germline_species.update(species_germlines)
    assert set(all_species) == germline_species
    assert set(all_species) == {"alpaca", "cow", "human", "mouse", "pig", "rabbit", "rat", "rhesus"}
