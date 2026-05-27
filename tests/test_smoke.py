"""Smoke test — confirms the package imports. Delete once real tests exist."""

import gem


def test_package_imports() -> None:
    assert gem.__version__ == "0.0.0"
