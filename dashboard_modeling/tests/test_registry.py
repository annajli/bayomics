from config.network_registry import NETWORKS, LEVEL_ORDER


def test_all_levels_present_and_ordered():
    assert LEVEL_ORDER == ["L1a", "L1b", "L2", "L3a", "L3b", "L4", "L5", "L6", "L_all"]
    assert set(NETWORKS) == set(LEVEL_ORDER)


def test_every_registry_file_exists():
    missing = []
    for level, meta in NETWORKS.items():
        for key in ("edges", "mapping", "png", "cv_edges"):
            if not meta[key].exists():
                missing.append(f"{level}.{key} -> {meta[key]}")
    assert not missing, "Missing artifacts:\n" + "\n".join(missing)
