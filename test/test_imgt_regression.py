"""IMGT germline reframe regression tests.

The corrected ALL.hmm/germlines in this fork remove spurious gap columns in
rhesus/mouse/rat germlines that used to shift the conserved IMGT anchors. These
tests assert the anchors land on their canonical IMGT columns: Cys23, Trp41,
Cys104 (and the documented FW1 frame for kappa).
"""
from anarci import number

SEQUENCES = {
    'heavy': 'QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYELVISCLDYWGQGTTLTVSS',
    'light': 'DVVMTQSPLSLPVTLGQPASISCRSSQSLVYSDGNTYLNWFQQRPGQSPRRLIYKVSNRDSGVPDRFSGSGSGTDFTLKISRVEAEDVGVYYCMQGTFGQGTKVEIK',
}


def test_imgt_conserved_anchors():
    for name, sequence in SEQUENCES.items():
        numbering, _chain_type = number(sequence, 'imgt')
        numbered = {key: aa for key, aa in numbering}
        assert numbered.get((23, ' ')) == 'C', name
        assert numbered.get((41, ' ')) == 'W', name
        assert numbered.get((104, ' ')) == 'C', name
        reproduced = ''.join(aa for _, aa in numbering if aa != '-')
        assert reproduced == sequence, name


def test_imgt_light_fw1_frame():
    """Human kappa FW1 must read 21:I 22:S 23:C 24:R.

    With the misframed germlines it was 21:- 22:I 23:S 24:C (GitHub issues
    #17, #24 pattern).
    """
    numbering, chain_type = number(SEQUENCES['light'], 'imgt')
    assert chain_type == 'L'
    by_position = {position: aa for (position, _insertion), aa in numbering}
    assert ''.join(by_position[position] for position in range(21, 25)) == 'ISCR'
