"""Tests for DEFAULT_DB_PATH resolution: the installer-managed location in
seamm.ini's [thermochemistry] section if set, else the bundled package path.
Uses monkeypatch to point at a throwaway ini file rather than touching the
real ~/.seamm.d/seamm.ini.
"""

from pathlib import Path

from seamm_thermochemistry import db as db_module


def _write_ini(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def test_falls_back_to_bundled_path_when_ini_missing(tmp_path, monkeypatch):
    monkeypatch.setattr(db_module, "_SEAMM_INI_PATH", tmp_path / "does_not_exist.ini")
    assert db_module._resolve_default_db_path() == db_module._BUNDLED_DB_PATH


def test_falls_back_when_no_thermochemistry_section(tmp_path, monkeypatch):
    ini = tmp_path / "seamm.ini"
    _write_ini(ini, "[SEAMM]\nroot = ~/SEAMM\n")
    monkeypatch.setattr(db_module, "_SEAMM_INI_PATH", ini)
    assert db_module._resolve_default_db_path() == db_module._BUNDLED_DB_PATH


def test_falls_back_when_database_path_key_missing(tmp_path, monkeypatch):
    ini = tmp_path / "seamm.ini"
    _write_ini(ini, "[thermochemistry]\nsome-other-key = value\n")
    monkeypatch.setattr(db_module, "_SEAMM_INI_PATH", ini)
    assert db_module._resolve_default_db_path() == db_module._BUNDLED_DB_PATH


def test_uses_installer_path_when_set(tmp_path, monkeypatch):
    installed = (
        tmp_path / "SEAMM" / "Parameters" / "thermochemistry" / "thermochemistry.db"
    )
    ini = tmp_path / "seamm.ini"
    _write_ini(ini, f"[thermochemistry]\ndatabase-path = {installed}\n")
    monkeypatch.setattr(db_module, "_SEAMM_INI_PATH", ini)
    assert db_module._resolve_default_db_path() == installed


def test_expands_user_in_database_path(tmp_path, monkeypatch):
    ini = tmp_path / "seamm.ini"
    _write_ini(
        ini, "[thermochemistry]\ndatabase-path = ~/some/path/thermochemistry.db\n"
    )
    monkeypatch.setattr(db_module, "_SEAMM_INI_PATH", ini)
    result = db_module._resolve_default_db_path()
    assert result == Path("~/some/path/thermochemistry.db").expanduser().resolve()


def test_tolerates_percent_signs_elsewhere_in_the_file(tmp_path, monkeypatch):
    # seamm.ini can contain literal "%" in unrelated values/comments (e.g.
    # ORCA's %pal, "memory-factor = 90%"); interpolation=None must mean
    # these never raise, even though we only ever read [thermochemistry].
    target = tmp_path / "thermochemistry.db"
    ini = tmp_path / "seamm.ini"
    _write_ini(
        ini,
        f"[orca-step]\nmemory-factor = 90%\n\n"
        f"[thermochemistry]\ndatabase-path = {target}\n",
    )
    monkeypatch.setattr(db_module, "_SEAMM_INI_PATH", ini)
    assert db_module._resolve_default_db_path() == target.resolve()
