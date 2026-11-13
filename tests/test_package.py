import ecoli_fba


def test_repo_paths_point_at_this_repository():
    assert (ecoli_fba.REPO_ROOT / "pyproject.toml").is_file()
    assert ecoli_fba.DATA_DIR == ecoli_fba.REPO_ROOT / "data"
    assert ecoli_fba.RESULTS_DIR == ecoli_fba.REPO_ROOT / "results"
