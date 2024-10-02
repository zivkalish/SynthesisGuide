import drawer, builder
from rdkit.Chem import MolFromSmiles, MolToSmiles
from reactions.halogenation import halogenation
from reactions.E2 import e2_reaction
from reactions.ozonolysis import ozonolysis


def main():
    s = "CCCC"
    mol = MolFromSmiles(s)
    # print(MolToSmiles(mol, kekuleSmiles=True, canonical=False))
    molecule = builder.molecule2graph(mol)
    drawer.draw_graph_as_mol(molecule)
    molecule = halogenation(molecule, "br")
    drawer.draw_graph_as_mol(molecule)
    molecule = e2_reaction(molecule, zaitsev=True)
    drawer.draw_graph_as_mol(molecule)
    molecule = ozonolysis(molecule, oxidative=True)[0]
    drawer.draw_graph_as_mol(molecule)


if __name__ == '__main__':
    main()
