# -*- coding: utf-8 -*-
"""In a shared installation the installer never downloads the database."""

import sys

import pytest

pytest.importorskip("seamm_installer")

from seamm_thermochemistry import installer as installer_module  # noqa: E402


@pytest.fixture
def inst(tmp_path, monkeypatch):
    monkeypatch.setattr(sys, "argv", ["seamm-thermochemistry-installer", "show"])
    monkeypatch.setattr(
        installer_module.Installer,
        "shared_codes",
        property(lambda self: True),
        raising=False,  # older seamm-manager has no shared_codes
    )
    installer = installer_module.Installer()
    db = tmp_path / "shared.db"
    monkeypatch.setattr(installer, "_configured_database_path", lambda: db)
    return installer, db


def test_shared_install_and_update_never_download(inst, monkeypatch, capsys):
    installer, db = inst

    def boom(*args, **kwargs):
        raise AssertionError("downloaded in a shared installation")

    monkeypatch.setattr(installer_module, "Zenodo", boom)
    installer.install()
    assert "Install it there first" in capsys.readouterr().out
    db.write_bytes(b"")
    installer.update()
    assert "uses the shared reference database" in capsys.readouterr().out
    installer.uninstall()
    assert db.exists()  # not deleted
