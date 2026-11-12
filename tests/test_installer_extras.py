# -*- coding: utf-8 -*-
"""The installer's official fragment store, official reference library and ATcT
(P5d, R3): the store is installed
from its Zenodo concept, verified, replaced whole by a newer published version
and kept when edited; ATcT's newest version is fetched when not cached. Zenodo
and ATcT are never contacted: the records and downloads are stand-ins.

The installer object is built without seamm_installer's InstallerBase.__init__
(which needs conda on the PATH, so the other installer tests skip in CI)."""

import hashlib
from pathlib import Path

import pytest

pytest.importorskip("seamm_installer")

from seamm_thermochemistry import atct as A  # noqa: E402
from seamm_thermochemistry import installer as installer_module  # noqa: E402
from seamm_thermochemistry import FragmentStore  # noqa: E402


def md5(data):
    return hashlib.md5(data).hexdigest()


class Configuration:
    """A stand-in for seamm_installer's configuration (seamm.ini)."""

    def __init__(self):
        self.sections = {}

    def section_exists(self, section):
        return section in self.sections

    def add_section(self, section):
        self.sections.setdefault(section, {})

    def get_values(self, section):
        return dict(self.sections.get(section, {}))

    def set_value(self, section, key, value):
        self.sections.setdefault(section, {})[key] = value

    def save(self):
        pass


def official_store(path, version):
    """A real (tiny) official store file, returned as bytes."""
    with FragmentStore(path, create=True, layer="official") as s:
        s.set_meta(version=version)
    return Path(path).read_bytes()


class Record(dict):
    def __init__(self, content, filename="fragments-official.db", checksum=None):
        super().__init__(
            id=24139306,
            files=[{"key": filename, "checksum": "md5:" + (checksum or md5(content))}],
        )
        self.content = content

    def download_file(self, filename, path):
        Path(path).write_bytes(self.content)


class Installer(installer_module.Installer):
    """The installer with a plain attribute for the configuration (the base
    class's is a read-only property over seamm.ini)."""

    configuration = None


@pytest.fixture()
def inst(tmp_path, monkeypatch):
    installer = Installer.__new__(Installer)
    installer.section = "thermochemistry"
    installer.configuration = Configuration()
    installer.options = None
    params = tmp_path / "SEAMM" / "Parameters" / "thermochemistry"
    db = params / "thermochemistry.db"
    monkeypatch.setattr(installer, "_configured_database_path", lambda: db)
    v1 = official_store(tmp_path / "v1.db", "2026.11.11")
    v2 = official_store(tmp_path / "v2.db", "2026.12.1")
    published = {md5(v1), md5(v2)}
    monkeypatch.setattr(installer, "_published_file_md5s", lambda r, f: published)
    installer._versions = {"v1": v1, "v2": v2}
    return installer, params


def test_install_then_up_to_date(inst, monkeypatch, capsys):
    installer, params = inst
    monkeypatch.setattr(
        installer, "_official_record", lambda: Record(installer._versions["v1"])
    )
    assert installer.install_official_store() == "installed"
    path = params / "fragments-official.db"
    assert path.read_bytes() == installer._versions["v1"]
    assert installer._official_version(path) == "2026.11.11"
    assert installer.configuration.get_values("thermochemistry")[
        "fragment-store-md5"
    ] == md5(installer._versions["v1"])
    assert installer.install_official_store() == "up to date"
    assert "version 2026.11.11" in capsys.readouterr().out


def test_a_newer_version_replaces_the_old_with_a_backup(inst, monkeypatch):
    installer, params = inst
    monkeypatch.setattr(
        installer, "_official_record", lambda: Record(installer._versions["v1"])
    )
    installer.install_official_store()
    monkeypatch.setattr(
        installer, "_official_record", lambda: Record(installer._versions["v2"])
    )
    assert installer.install_official_store() == "updated"
    path = params / "fragments-official.db"
    assert installer._official_version(path) == "2026.12.1"
    assert len(list(params.glob("fragments-official.db.bak-*"))) == 1


def test_an_edited_store_is_kept_unless_forced(inst, monkeypatch, capsys):
    installer, params = inst
    path = params / "fragments-official.db"
    params.mkdir(parents=True)
    path.write_bytes(b"edited by hand")
    monkeypatch.setattr(
        installer, "_official_record", lambda: Record(installer._versions["v2"])
    )
    assert installer.install_official_store() == "kept"
    assert path.read_bytes() == b"edited by hand"
    assert "not a published version" in capsys.readouterr().out
    assert installer.install_official_store(force=True) == "updated"
    assert installer._official_version(path) == "2026.12.1"


def test_a_bad_download_never_replaces_the_store(inst, monkeypatch, capsys):
    installer, params = inst
    monkeypatch.setattr(
        installer, "_official_record", lambda: Record(installer._versions["v1"])
    )
    installer.install_official_store()
    bad = Record(installer._versions["v2"], checksum="0" * 32)
    monkeypatch.setattr(installer, "_official_record", lambda: bad)
    assert installer.install_official_store() == "failed"
    path = params / "fragments-official.db"
    assert installer._official_version(path) == "2026.11.11"
    assert not path.with_name(path.name + ".download").exists()
    assert "Checksum mismatch" in capsys.readouterr().out


def test_no_network_is_reported_not_raised(inst, monkeypatch, capsys):
    installer, _ = inst

    def offline():
        raise ConnectionError("no route to zenodo.org")

    monkeypatch.setattr(installer, "_official_record", offline)
    assert installer.install_official_store() == "failed"
    assert "no route" in capsys.readouterr().out


def reference_library(path, version):
    """A real (tiny) official reference library file, returned as bytes."""
    pytest.importorskip("rdkit")
    from seamm_thermochemistry.references import ReferenceLibrary

    library = ReferenceLibrary(path, create=True, layer="official")
    library.set_meta("version", version)
    library.close()
    return Path(path).read_bytes()


def test_reference_library_install_update_and_kept(inst, tmp_path, monkeypatch):
    """The official reference library is installed beside the ATcT directory,
    where the Thermochemistry step reads it, replaced by a newer published
    version with a backup, and kept when it is not a published version."""
    installer, params = inst
    v1 = reference_library(tmp_path / "r1.db", "2026.11.12")
    v2 = reference_library(tmp_path / "r2.db", "2026.12.1")
    monkeypatch.setattr(
        installer, "_published_file_md5s", lambda r, f: {md5(v1), md5(v2)}
    )
    name = "references-official.db"
    monkeypatch.setattr(installer, "_references_record", lambda: Record(v1, name))
    assert installer.install_reference_library() == "installed"
    path = params / name
    assert path.parent == installer._atct_base().parent
    assert installer._references_version(path) == "2026.11.12"
    assert installer.configuration.get_values("thermochemistry")[
        "reference-library-md5"
    ] == md5(v1)
    assert installer.install_reference_library() == "up to date"
    monkeypatch.setattr(installer, "_references_record", lambda: Record(v2, name))
    assert installer.install_reference_library() == "updated"
    assert installer._references_version(path) == "2026.12.1"
    assert len(list(params.glob(name + ".bak-*"))) == 1
    path.write_bytes(b"edited by hand")
    assert installer.install_reference_library() == "kept"
    assert installer.install_reference_library(force=True) == "updated"


FIXTURE = Path(__file__).parent / "data" / "atct_fixture.json"


def test_install_atct(inst, monkeypatch, capsys):
    installer, params = inst
    real = A.download_all

    def from_fixture(base, **kwargs):
        return real(base, source=FIXTURE, **kwargs)

    monkeypatch.setattr(A, "download_all", from_fixture)
    version = installer.install_atct()
    assert version and (params / "ATcT" / version / "atct_all.json").exists()
    assert "Downloaded ATcT" in capsys.readouterr().out
    # Again: a version is cached, so nothing is downloaded, and newer ones are
    # added only on purpose
    calls = []
    monkeypatch.setattr(A, "download_all", lambda *a, **k: calls.append(1))
    assert installer.install_atct() == version and calls == []
    assert "only on purpose" in capsys.readouterr().out
    # On purpose, with the newest version already there: not an error
    monkeypatch.setattr(A, "download_all", from_fixture)
    assert installer.install_atct(newest=True) == version
    assert "already in" in capsys.readouterr().out


def test_install_atct_refused(inst, monkeypatch, capsys):
    installer, _ = inst

    def refused(base, **kwargs):
        raise A.ATcTError("ATcT refused the request: HTTP 403.")

    monkeypatch.setattr(A, "download_all", refused)
    assert installer.install_atct() is None
    out = capsys.readouterr().out
    assert "HTTP 403" in out and "download --from FILE" in out


# -- the ATcT cache ---------------------------------------------------------------


def test_download_all_exist_ok(tmp_path):
    directory = A.download_all(tmp_path, source=FIXTURE)
    with pytest.raises(A.ATcTError, match="already"):
        A.download_all(tmp_path, source=FIXTURE)
    assert A.download_all(tmp_path, source=FIXTURE, exist_ok=True) == directory
    assert not list(directory.glob("atct_all.json.download*"))


def test_covariance_cache_merges_on_save(tmp_path):
    """Two processes (here: two ATcT objects) caching different blocks keep both."""
    A.download_all(tmp_path, source=FIXTURE)
    first = A.ATcT(base=tmp_path)
    second = A.ATcT(base=tmp_path)
    first.covariance_cache["blocks"].append({"ids": ["a", "b"], "matrix": [[1]]})
    first._save_covariance()
    second.covariance_cache["blocks"].append({"ids": ["c"], "matrix": [[2]]})
    second._save_covariance()
    third = A.ATcT(base=tmp_path)
    assert {tuple(b["ids"]) for b in third.covariance_cache["blocks"]} == {
        ("a", "b"),
        ("c",),
    }
