from __future__ import annotations

from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors


@dataclass(frozen=True)
class RuleResult:
    """Outcome and summary explanation for a drug-likeness rule."""

    name: str
    passed: bool
    details: str
    violations: int = 0


def lipinski(mw: float, logp: float, hbd: int, hba: int) -> RuleResult:
    """Evaluate Lipinski's rule of five using common threshold checks."""
    failures = [
        label
        for label, bad in (
            ("MW > 500", mw > 500),
            ("cLogP > 5", logp > 5),
            ("HBD > 5", hbd > 5),
            ("HBA > 10", hba > 10),
        )
        if bad
    ]
    passed = len(failures) <= 1
    return RuleResult(
        "Lipinski",
        passed,
        "No violations" if not failures else "; ".join(failures),
        len(failures),
    )


def veber(rotatable_bonds: int, tpsa: float) -> RuleResult:
    """Evaluate Veber's oral availability criteria."""
    failures = [
        label
        for label, bad in (
            ("rotatable bonds > 10", rotatable_bonds > 10),
            ("TPSA > 140 Å²", tpsa > 140),
        )
        if bad
    ]
    return RuleResult(
        "Veber",
        not failures,
        "No violations" if not failures else "; ".join(failures),
        len(failures),
    )


def ghose(mw: float, logp: float, molar_refractivity: float, atoms: int) -> RuleResult:
    """Evaluate the Ghose filter.

    The Ghose filter is commonly stated as: 160 <= MW <= 480, -0.4 <= LogP <= 5.6,
    40 <= MR <= 130, and 20 <= atom count <= 70.
    """
    failures = [
        label
        for label, bad in (
            ("MW outside 160–480", not 160 <= mw <= 480),
            ("cLogP outside -0.4–5.6", not -0.4 <= logp <= 5.6),
            ("molar refractivity outside 40–130", not 40 <= molar_refractivity <= 130),
            ("atom count outside 20–70", not 20 <= atoms <= 70),
        )
        if bad
    ]
    return RuleResult(
        "Ghose",
        not failures,
        "All ranges satisfied" if not failures else "; ".join(failures),
        len(failures),
    )


def egan(logp: float, tpsa: float) -> RuleResult:
    """Evaluate the Egan filter."""
    failures = [
        label
        for label, bad in (
            ("cLogP > 5.88", logp > 5.88),
            ("TPSA > 131.6 Å²", tpsa > 131.6),
        )
        if bad
    ]
    return RuleResult(
        "Egan",
        not failures,
        "No violations" if not failures else "; ".join(failures),
        len(failures),
    )


def muegge(mw: float, logp: float, tpsa: float, hbd: int, hba: int, rotors: int, rings: int) -> RuleResult:
    """Evaluate a practical Muegge-type drug-likeness filter."""
    failures = [
        label
        for label, bad in (
            ("MW outside 200–600", not 200 <= mw <= 600),
            ("cLogP outside -2–5", not -2 <= logp <= 5),
            ("TPSA > 150 Å²", tpsa > 150),
            ("HBD > 5", hbd > 5),
            ("HBA > 10", hba > 10),
            ("rotatable bonds > 15", rotors > 15),
            ("rings > 7", rings > 7),
        )
        if bad
    ]
    return RuleResult(
        "Muegge",
        not failures,
        "All ranges satisfied" if not failures else "; ".join(failures),
        len(failures),
    )


def evaluate_rules(molecule: Chem.Mol) -> list[RuleResult]:
    """Compute the main oral drug-likeness filters for a valid molecule."""
    if molecule is None:
        raise ValueError("A valid RDKit molecule is required.")

    mw = float(Descriptors.MolWt(molecule))
    logp = float(Descriptors.MolLogP(molecule))
    tpsa = float(rdMolDescriptors.CalcTPSA(molecule))
    hbd = int(Lipinski.NumHDonors(molecule))
    hba = int(Lipinski.NumHAcceptors(molecule))
    rotors = int(Lipinski.NumRotatableBonds(molecule))
    rings = int(rdMolDescriptors.CalcNumRings(molecule))
    atoms = int(molecule.GetNumAtoms())
    refractivity = float(Descriptors.MolMR(molecule))

    return [
        lipinski(mw, logp, hbd, hba),
        veber(rotors, tpsa),
        ghose(mw, logp, refractivity, atoms),
        egan(logp, tpsa),
        muegge(mw, logp, tpsa, hbd, hba, rotors, rings),
    ]
