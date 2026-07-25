"""Derive DfH(X, 0K) from the 298 K experimental data already in the store.

The master reference sheet only has ΔfH°(0K) filled in for H and He; every
other element has ΔfH°(298K) plus the two H(298K)-H(0K) corrections needed
to convert it, via the standard thermochemical identity:

    DfH(X, 0K) = DfH(X, 298K) - [H298-H0]_atom + [H298-H0]_standard_state

`h298_minus_h0_std_state` is tabulated *per formula unit of the standard
state*, not per atom -- for a diatomic standard state (H2, N2, O2, F2, Cl2,
Br2) it must be divided by 2 to get a per-atom contribution; for the
monatomic/solid standard states used by every other element in range, no
division is needed. Verified against the two elements the master sheet
already gives directly: H (218 - 6.197 + 8.468/2 = 216.03, matches
216.034) and He (0 - 6.197 + 6.197 = 0, matches 0), both to within the
sheet's own rounding.
"""

__all__ = ["derive_dfH0_0K"]

# Elements whose standard state is a diatomic gas (X2) -- everything else in
# the periodic table's H-Kr range has a monatomic formula unit (solid metal,
# noble gas, C(s,gr), etc.), so no factor is needed.
_DIATOMIC_STANDARD_STATE = {"H", "N", "O", "F", "Cl", "Br"}


def derive_dfH0_0K(db, *, overwrite=False):
    """Fill in `element.dfH0_0K` for every element where it's derivable.

    Parameters
    ----------
    db : ThermoDB
    overwrite : bool
        If False (default), skip elements that already have `dfH0_0K` --
        the master sheet's own H/He values are treated as authoritative
        over a re-derivation (though the derivation reproduces them to
        within rounding; see module docstring).

    Returns
    -------
    list of str
        Symbols of the elements filled in.
    """
    filled = []
    for el in db.elements():
        if not overwrite and el["dfH0_0K"] is not None:
            continue
        dfH298 = el["dfH0_298K"]
        h298h0_atom = el["h298_minus_h0_atom"]
        h298h0_std = el["h298_minus_h0_std_state"]
        if dfH298 is None or h298h0_atom is None or h298h0_std is None:
            continue

        divisor = 2 if el["symbol"] in _DIATOMIC_STANDARD_STATE else 1
        dfH0_0K = dfH298 - h298h0_atom + h298h0_std / divisor

        db.add_element(
            el["atomic_number"],
            el["symbol"],
            multiplicity=el["multiplicity"],
            term_symbol=el["term_symbol"],
            standard_state=el["standard_state"],
            dfH0_0K=dfH0_0K,
            dfH0_298K=el["dfH0_298K"],
            dfH0_298K_stderr=el["dfH0_298K_stderr"],
            h298_minus_h0_atom=el["h298_minus_h0_atom"],
            h298_minus_h0_std_state=el["h298_minus_h0_std_state"],
            s298_gas=el["s298_gas"],
            s298_gas_stderr=el["s298_gas_stderr"],
            reference=el["reference"],
            reference_note=(
                f"{el['reference_note']} (0K derived from 298K)"
                if el["reference_note"]
                else "0K derived from 298K"
            ),
        )
        filled.append(el["symbol"])
    return filled
