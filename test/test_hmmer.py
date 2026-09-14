"""Tests for the pyhmmer-backed alignment step (run_pyhmmer)."""

import pytest

from anarci import run_pyhmmer

SEQUENCE = "QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYELVISCLDYWGQGTTLTVSS"

# Alignment against QVQLQQSGA-ELARPGASVKMSCKASGYTF----TRYTMHWVKQRPGQGLEWIGYINPS--
#                  RGYTNYNQKFK-DKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYELVISCLDYWGQGTTLTVSS
STATE_TYPES = "mmmmmmmmmdmmmmmmmmmmmmmmmmmmmmddddmmmmmmmmmmmmmmmmmmmmmmmmmddmmmmmmmmmmmdmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmiiiimmmmmmmmmmmmmmmmmm"


def test_run_pyhmmer():
    hit_table, state_vectors, top_descriptions = run_pyhmmer(
        [("id", SEQUENCE)],
        hmm_database="ALL",
        ncpu=None,
        bit_score_threshold=80,
        hmmer_species=["human", "mouse"],
    )[0]

    state_vector = state_vectors[0]

    hmm_states = [hmm_state for (hmm_state, _state_type), _sequence_index in state_vector]
    state_types = "".join(state_type for (_hmm_state, state_type), _sequence_index in state_vector)
    sequence_reproduced = "".join(SEQUENCE[sequence_index] if sequence_index is not None else "" for _key, sequence_index in state_vector)

    assert sequence_reproduced == SEQUENCE
    assert state_types == STATE_TYPES
    assert hmm_states == list(range(1, 111)) + [111, 111, 111, 111, 111] + list(range(112, 129))

    assert len(hit_table) == 3
    assert hit_table[0] == ["id", "description", "evalue", "bitscore", "bias", "query_start", "query_end"]
    assert hit_table[1][0] == "mouse_H"
    assert hit_table[1][1] == ""
    # Values reflect the corrected ALL.hmm (IMGT germline reframe + pyhmmer 0.12).
    assert hit_table[1][2] == pytest.approx(2.7e-60)
    assert hit_table[1][3] == pytest.approx(193.1, rel=1e-2)
    assert hit_table[1][4] == pytest.approx(3.59, rel=1e-2)
    assert hit_table[1][5] == 0
    assert hit_table[1][6] == 124

    assert len(top_descriptions) == 1
    top = top_descriptions[0]
    assert top["id"] == "mouse_H"
    assert top["description"] == ""
    assert top["evalue"] == pytest.approx(2.7e-60)
    assert top["bitscore"] == pytest.approx(193.1, rel=1e-2)
    assert top["bias"] == pytest.approx(3.59, rel=1e-2)
    assert top["query_start"] == 0
    assert top["query_end"] == 124
    assert top["chain_type"] == "H"
    assert top["species"] == "mouse"
