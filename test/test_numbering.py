"""IMGT numbering tests for the single-sequence number() helper."""

from anarci import number

HEAVY_SEQUENCE = "QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYELVISCLDYWGQGTTLTVSS"
LIGHT_SEQUENCE = "DVVMTQSPLSLPVTLGQPASISCRSSQSLVYSDGNTYLNWFQQRPGQSPRRLIYKVSNRDSGVPDRFSGSGSGTDFTLKISRVEAEDVGVYYCMQGTFGQGTKVEIK"
# Heavy chain with a few mismatching residues at the N and C termini.
HEAVY_INCOMPLETE_SEQUENCE = (
    "AAQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYELVISCLDYWGQGTTLTVAA"
)

HEAVY_ALIGNED = "QVQLQQSGA-ELARPGASVKMSCKASGYTF----TRYTMHWVKQRPGQGLEWIGYINPS--RGYTNYNQKFK-DKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYELVISCLDYWGQGTTLTVSS"
LIGHT_ALIGNED = "DVVMTQSPLSLPVTLGQPASISCRSSQSLVYS-DGNTYLNWFQQRPGQSPRRLIYKV-------SNRDSGVP-DRFSGSG--SGTDFTLKISRVEAEDVGVYYCMQ---------GTFGQGTKVEIK"
HEAVY_INCOMPLETE_ALIGNED = (
    "AAQLQQSGA-ELARPGASVKMSCKASGYTF----TRYTMHWVKQRPGQGLEWIGYINPS--RGYTNYNQKFK-DKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYELVISCLDYWGQGTTLTVAA"
)


def _positions(numbering):
    return [f"{position}{insertion}".strip() for (position, insertion), _amino_acid in numbering]


def _aligned(numbering):
    return "".join(amino_acid for _key, amino_acid in numbering)


def _reproduced(numbering):
    return "".join(amino_acid for _key, amino_acid in numbering if amino_acid != "-")


def test_imgt_heavy():
    numbering, chain_type = number(HEAVY_SEQUENCE, "imgt")

    assert chain_type == "H"
    expected_positions = [str(i) for i in range(1, 111)] + ["111", "111A", "111B", "112B", "112A"] + [str(i) for i in range(112, 129)]
    assert _positions(numbering) == expected_positions
    assert _aligned(numbering) == HEAVY_ALIGNED
    assert _reproduced(numbering) == HEAVY_SEQUENCE


def test_imgt_light():
    numbering, chain_type = number(LIGHT_SEQUENCE, "imgt")

    assert chain_type == "L"
    assert _positions(numbering) == [str(i) for i in range(1, 128)]
    assert _aligned(numbering) == LIGHT_ALIGNED
    assert _reproduced(numbering) == LIGHT_SEQUENCE


def test_imgt_heavy_incomplete():
    """A few mismatching residues at the N and C termini stay part of the aligned domain."""
    numbering, _chain_type = number(HEAVY_INCOMPLETE_SEQUENCE, "imgt")

    assert _aligned(numbering) == HEAVY_INCOMPLETE_ALIGNED
    assert _reproduced(numbering) == HEAVY_INCOMPLETE_SEQUENCE
