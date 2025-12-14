from __future__ import annotations

from sc_browser.core.filter_state import FilterState


def test_filter_state_to_from_dict_roundtrip():
    st = FilterState(
        dataset_name="ds1",
        view_id="expression",
        genes=["G1", "G2"],
        clusters=["A"],
        conditions=["C1"],
        samples=["S1"],
        cell_types=["T"],
        embedding="X_umap",
        merge_genes=True,
        split_by_condition=True,
        is_3d=True,
        dim_x=1,
        dim_y=0,
        dim_z=2,
        color_scale="plasma",
    )

    raw = st.to_dict()
    rebuilt = FilterState.from_dict(raw)

    assert rebuilt == st
