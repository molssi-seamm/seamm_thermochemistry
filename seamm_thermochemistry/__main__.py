# !/usr/bin/env python
# -*- coding: utf-8 -*-

"""Handle the installation of the seamm_thermochemistry reference database."""

from .installer import Installer


def run():
    """Handle the extra installation needed: download/update/remove the
    thermochemistry reference database and register its location in
    seamm.ini."""

    installer = Installer()
    installer.run()


if __name__ == "__main__":
    run()
