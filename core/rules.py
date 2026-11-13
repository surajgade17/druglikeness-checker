"""Drug-likeness rule implementations."""
from __future__ import annotations

from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors


@dataclass(frozen=True)
class RuleResult:
    """Outcome and explanation for one filter."""
    name: str
    passed: bool
    details: str
    violations: int = 0


def _result(name: str, failures: list[str], allow_one: bool = False) -> RuleResult:
    return RuleResult(name, len(failures) <= (1 if allow_one else 0), "No violations" if not failures else "; ".join(failures), len(failures))


def lipinski(mw: float, logp: float, hbd: int, hba: int) -> RuleResult:
    """Lipinski rule of five; at most one violation is allowed."""
    failures = [label for label, bad in (("MW > 500", mw > 500), ("cLogP > 5", logp > 5), ("HBD > 5", hbd > 5), ("HBA > 10", hba > 10)) if bad]
    return _result("Lipinski", failures, allow_one=True)


def veber(rotatable_bonds: int, tpsa: float) -> RuleResult:
    """Veber oral bioavailability criteria."""
    failures = [label for label, bad in (("rotatable bonds > 10", rotatable_bonds > 10), ("TPSA > 140 Å²", tpsa > 140)) if bad]
    return _result("Veber", failures)


def ghose(mw: float, logp: float, molar_refractivity: float, atoms: int) -> RuleResult:
    """Ghose filter ranges."""
    failures = [label for label, bad in (("MW outside 160–480", not 160 <= mw <= 480), ("cLogP outside -0.4–5.6", not -0.4 <= logp <= 5.6), ("molar refractivity outside 40–130", not 40 <= molar_refractivity <= 130), ("atom count outside 20–70", not 20 <= atoms <= 70)) if bad]
    return _result("Ghose", failures)


def egan(logp: float, tpsa: float) -> RuleResult:
    """Egan filter."""
    failures = [label for label, bad in (("cLogP > 5.88", logp > 5.88), ("TPSA > 131.6 Å²", tpsa > 131.6)) if bad]
    return _result("Egan", failures)


def muegge(mw: float, logp: float, tpsa: float, hbd: int, hba: int, rotors: int, rings: int) -> RuleResult:
    """Commonly used Muegge filter bounds."""
    failures = [label for label, bad in (("MW outside 200–600", not 200 <= mw <= 600), ("cLogP outside -2–5", not -2 <= logp <= 5), ("TPSA > 150 Å²", tpsa > 150), ("HBD > 5", hbd > 5), ("HBA > 10", hba > 10), ("rotatable bonds > 15", rotors > 15), ("rings > 7", rings > 7)) if bad]
    return _result("Muegge", failures)


def evaluate_rules(molecule: Chem.Mol) -> list[RuleResult]:
    """Evaluate five drug-likeness filters for a valid molecule."""
    if molecule is None:
        raise ValueError("A valid RDKit molecule is required.")
    mw = float(Descriptors.MolWt(molecule))
    logp = float(Descriptors.MolLogP(molecule))
    tpsa = float(rdMolDescriptors.CalcTPSA(molecule))
    hbd, hba = int(Lipinski.NumHDonors(molecule)), int(Lipinski.NumHAcceptors(molecule))
    rotors = int(Lipinski.NumRotatableBonds(molecule))
    rings = int(rdMolDescriptors.CalcNumRings(molecule))
    atoms, refractivity = int(molecule.GetNumAtoms()), float(Descriptors.MolMR(molecule))
    return [lipinski(mw, logp, hbd, hba), veber(rotors, tpsa), ghose(mw, logp, refractivity, atoms), egan(logp, tpsa), muegge(mw, logp, tpsa, hbd, hba, rotors, rings)]
