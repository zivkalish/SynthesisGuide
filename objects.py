from dataclasses import dataclass
from rdkit import Chem

@dataclass(frozen=True)
class Atom:
    idx: int
    atomic_num: int
    symbol: str
    charge: int

    @classmethod
    def from_chem(cls, atom: Chem.Atom):
        return cls(
            idx=atom.GetIdx(),
            atomic_num=atom.GetAtomicNum(),
            symbol=atom.GetSymbol(),
            charge=atom.GetFormalCharge(),
        )

    def to_chem(self) -> Chem.Atom:
        atom = Chem.Atom(self.atomic_num)
        atom.SetFormalCharge(self.charge)
        return atom
