# -*- coding: utf-8 -*-
"""A fragment store from newer code (a higher ``PRAGMA user_version``) is never
written by this code, and reading it is warned about once: older code would not
maintain what a newer schema adds (e.g. a quarantine column)."""

import logging
import sqlite3

import pytest

from seamm_thermochemistry import FragmentRecord, FragmentStore
from seamm_thermochemistry import fragment_store as F
from seamm_thermochemistry.geometry_store import SCHEMA_VERSION


@pytest.fixture()
def newer(tmp_path):
    """A store with one record, then marked as one schema version newer."""
    path = tmp_path / "fragments.db"
    with FragmentStore(path, create=True) as store:
        store.put(
            FragmentRecord.make(
                "C", E=-100.0, ZPE=1.0, H=-90.0, model_chemistry="Test:X"
            )
        )
    db = sqlite3.connect(path)
    db.execute(f"PRAGMA user_version = {SCHEMA_VERSION + 1}")
    db.commit()
    db.close()
    F._warned_newer.clear()
    return path


def test_a_newer_store_is_not_written(newer):
    before = newer.read_bytes()
    with pytest.raises(F.NewerSchemaError, match="not written to"):
        FragmentStore(newer, create=True)
    # The tools that write through an opened store are refused too
    store = FragmentStore(newer)
    store._readonly = False
    with pytest.raises(F.NewerSchemaError):
        store._connect(create=True)
    assert newer.read_bytes() == before
    assert F.file_schema_version(newer) == SCHEMA_VERSION + 1


def test_a_newer_store_is_read_with_one_warning(newer, caplog):
    with caplog.at_level(logging.WARNING, logger=F.__name__):
        for _ in range(2):
            with FragmentStore(newer) as store:
                assert store.schema_version() == SCHEMA_VERSION + 1
                assert len(store.records("C")) == 1
    warnings = [r for r in caplog.records if "newer than" in r.getMessage()]
    assert len(warnings) == 1


def test_a_newer_store_is_never_writable(newer, monkeypatch, tmp_path):
    ini = tmp_path / "seamm.ini"
    ini.write_text("[thermochemistry]\nfragment-store-writable = yes\n")
    monkeypatch.setattr("seamm_thermochemistry.db._SEAMM_INI_PATH", ini)
    writable, why = F.store_is_writable(newer)
    assert not writable
    assert f"schema version {SCHEMA_VERSION + 1}" in why


def test_current_and_new_stores_are_unaffected(tmp_path, monkeypatch):
    ini = tmp_path / "seamm.ini"
    ini.write_text("[thermochemistry]\nfragment-store-writable = yes\n")
    monkeypatch.setattr("seamm_thermochemistry.db._SEAMM_INI_PATH", ini)
    path = tmp_path / "fragments.db"
    assert F.store_is_writable(path)[0]  # no file yet
    with FragmentStore(path, create=True) as store:
        assert store.put(
            FragmentRecord.make("C", E=-1.0, ZPE=1.0, H=0.0, model_chemistry="Test:X")
        )
    assert F.file_schema_version(path) == SCHEMA_VERSION
    assert F.store_is_writable(path)[0]
    assert F.file_schema_version(tmp_path / "missing.db") in (None, 0)
