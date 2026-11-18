# -*- coding: utf-8 -*-
"""Zenodo errors in the installer (2026-11-18): transient errors (5xx, 429, an
empty reply, a dropped connection) are tried again; anything still failing is
reported and the installer exits with an error, rather than printing "Could not
install ..." and exiting 0, which left installations on old data unnoticed."""

import json

import pytest

pytest.importorskip("seamm_installer")

from seamm_thermochemistry import installer as installer_module  # noqa: E402

from tests import test_installer_extras as extras  # noqa: E402

Record = extras.Record
inst = extras.inst  # the shared fixture


@pytest.fixture(autouse=True)
def no_waiting(monkeypatch):
    waits = []
    monkeypatch.setattr(installer_module.time, "sleep", waits.append)
    return waits


def failing(errors, then):
    """A function that raises each of `errors` in turn, then returns `then`."""
    calls = []

    def function(*args):
        calls.append(args)
        if len(calls) <= len(errors):
            raise errors[len(calls) - 1]
        return then

    function.calls = calls
    return function


@pytest.mark.parametrize(
    "error, transient",
    [
        (RuntimeError("Error in get_latest_public_record: code = 500"), True),
        (RuntimeError("Error in get_file: code = 503"), True),
        (RuntimeError("Error in download_file: code = 429"), True),
        (json.JSONDecodeError("Expecting value", "", 0), True),
        (ConnectionError("reset"), True),
        (RuntimeError("Error in get_latest_public_record: code = 404"), False),
        (RuntimeError("Checksum mismatch downloading x"), False),
        (KeyError("files"), False),
    ],
)
def test_is_transient(error, transient):
    assert installer_module.is_transient(error) is transient


def test_transient_errors_are_tried_again(inst, no_waiting, capsys):
    installer, _ = inst
    error = RuntimeError("Error in get_latest_public_record: code = 500")
    function = failing([error, error], "record")
    assert installer._zenodo("the record", function, 1) == "record"
    assert len(function.calls) == 3
    assert no_waiting == [5, 15]
    assert "trying again in 5 s" in capsys.readouterr().out


def test_other_errors_and_the_last_try_are_raised(inst, no_waiting, monkeypatch):
    installer, _ = inst
    missing = RuntimeError("Error in get_latest_public_record: code = 404")
    function = failing([missing], "record")
    with pytest.raises(RuntimeError, match="404"):
        installer._zenodo("the record", function)
    assert len(function.calls) == 1 and no_waiting == []

    monkeypatch.setenv("SEAMM_ZENODO_ATTEMPTS", "2")
    error = RuntimeError("code = 502")
    function = failing([error, error, error], "record")
    with pytest.raises(RuntimeError, match="502"):
        installer._zenodo("the record", function)
    assert len(function.calls) == 2


class FlakyZenodo:
    """Zenodo whose record requests fail with 500 a given number of times."""

    failures = 0
    records = {}

    def get_latest_public_record(self, concept):
        if FlakyZenodo.failures > 0:
            FlakyZenodo.failures -= 1
            raise RuntimeError("Error in get_latest_public_record: code = 500")
        if concept not in FlakyZenodo.records:
            raise RuntimeError("Error in get_latest_public_record: code = 404")
        return FlakyZenodo.records[concept]


def test_a_brief_outage_still_installs(inst, monkeypatch):
    installer, params = inst
    v1 = installer._versions["v1"]
    FlakyZenodo.failures = 2
    FlakyZenodo.records = {installer.fragment_store_concept_id: Record(v1)}
    monkeypatch.setattr(installer_module, "Zenodo", FlakyZenodo)
    assert installer.install_official_store() == "installed"
    assert (params / "fragments-official.db").read_bytes() == v1
    assert installer.failures == []


def test_update_exits_with_an_error_and_says_what_failed(inst, monkeypatch, capsys):
    """Zenodo down throughout: every part is tried, each failure is named, and
    the installer exits 1 instead of 0."""
    installer, params = inst
    FlakyZenodo.failures = 10**6
    FlakyZenodo.records = {}
    monkeypatch.setattr(installer_module, "Zenodo", FlakyZenodo)
    monkeypatch.setattr(installer, "_shared_report", lambda: False)
    monkeypatch.setattr(installer, "install_atct", lambda newest=False: None)
    with pytest.raises(SystemExit) as exit_:
        installer.update()
    assert exit_.value.code == 1
    out = capsys.readouterr().out
    assert installer.failures == [
        "the thermochemistry database",
        "the official fragment store",
        "the official reference library",
        "the computed reference values",
    ]
    assert "Could not install or update the thermochemistry database, the " in out
    assert "seamm-thermochemistry-installer update' again later" in out


def test_update_without_errors_exits_normally(inst, monkeypatch):
    installer, params = inst
    monkeypatch.setattr(installer, "_shared_report", lambda: False)
    monkeypatch.setattr(installer, "_update_database", lambda: None)
    for name in (
        "install_official_store",
        "install_reference_library",
        "install_computed_values",
    ):
        monkeypatch.setattr(installer, name, lambda force=False: "up to date")
    monkeypatch.setattr(installer, "install_atct", lambda newest=False: None)
    installer.update()  # no SystemExit
    assert installer.failures == []
