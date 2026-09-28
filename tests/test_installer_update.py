# -*- coding: utf-8 -*-
"""The installer never overwrites or deletes a database with local changes."""

import hashlib
from pathlib import Path
import sys

import pytest

pytest.importorskip("seamm_installer")

from seamm_thermochemistry import installer as installer_module  # noqa: E402

OLD, LATEST, CURATED = b"published v1", b"published v2", b"curated edits"


def md5(data):
    return hashlib.md5(data).hexdigest()


class _Record(dict):
    def __init__(self, content, checksum=None):
        super().__init__(
            id=2,
            files=[
                {
                    "key": "thermochemistry.db",
                    "checksum": "md5:" + (checksum or md5(content)),
                }
            ],
        )
        self.content = content

    def download_file(self, filename, path):
        Path(path).write_bytes(self.content)


@pytest.fixture
def inst(tmp_path, monkeypatch):
    home = tmp_path / "home"
    (home / ".seamm.d").mkdir(parents=True)
    monkeypatch.setenv("HOME", str(home))  # never the real ~/.seamm.d/seamm.ini
    monkeypatch.setenv("SEAMM_ROOT", str(home / "SEAMM"))
    monkeypatch.setattr(sys, "argv", ["seamm-thermochemistry-installer", "show"])
    monkeypatch.setattr(
        installer_module.Installer,
        "shared_codes",
        property(lambda self: False),
        raising=False,
    )
    installer = installer_module.Installer()
    db = home / "SEAMM" / "Parameters" / "thermochemistry" / "thermochemistry.db"
    monkeypatch.setattr(installer, "_configured_database_path", lambda: db)
    monkeypatch.setattr(installer, "_latest_record", lambda: _Record(LATEST))
    monkeypatch.setattr(
        installer, "_published_md5s", lambda record: {md5(OLD), md5(LATEST)}
    )
    return installer, db


def backups(db):
    return sorted(db.parent.glob(db.name + ".bak-*"))


def test_install_then_up_to_date(inst, capsys):
    installer, db = inst
    installer.install()
    assert db.read_bytes() == LATEST
    installer.update()
    assert "up to date" in capsys.readouterr().out
    assert backups(db) == []


def test_update_replaces_an_older_published_copy(inst):
    installer, db = inst
    db.parent.mkdir(parents=True)
    db.write_bytes(OLD)
    installer.update()
    assert db.read_bytes() == LATEST
    assert [b.read_bytes() for b in backups(db)] == [OLD]


def test_update_keeps_local_changes_unless_forced(inst, capsys):
    installer, db = inst
    db.parent.mkdir(parents=True)
    db.write_bytes(CURATED)
    installer.update()
    assert db.read_bytes() == CURATED
    assert "has local changes" in capsys.readouterr().out
    installer.install()  # install never overwrites either
    assert db.read_bytes() == CURATED

    installer.options = type("O", (), {"force": True})()
    installer.update()
    assert db.read_bytes() == LATEST
    assert [b.read_bytes() for b in backups(db)] == [CURATED]


def test_failed_download_leaves_the_database_alone(inst, monkeypatch):
    installer, db = inst
    db.parent.mkdir(parents=True)
    db.write_bytes(OLD)
    monkeypatch.setattr(
        installer, "_latest_record", lambda: _Record(b"corrupt", checksum=md5(LATEST))
    )
    with pytest.raises(RuntimeError, match="Checksum mismatch"):
        installer.update()
    assert db.read_bytes() == OLD
    assert not db.with_name(db.name + ".download").exists()


def test_uninstall_keeps_local_changes(inst):
    installer, db = inst
    installer.install()
    db.write_bytes(CURATED)
    installer.uninstall()
    assert db.exists()
    installer.install()  # re-registers; the curated file is kept
    db.write_bytes(LATEST)  # now an unmodified published copy
    installer.uninstall()
    assert not db.exists()
