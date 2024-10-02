import drawer, builder
from rdkit.Chem import MolFromSmiles, MolToSmiles
from reactions.halogenation import halogenation

def main():
    s = "CCCC"
    mol = MolFromSmiles(s)
    # print(MolToSmiles(mol, kekuleSmiles=True, canonical=False))
    molecule = builder.molecule2graph(mol)
    # drawer.draw_graph(graph)
    drawer.draw_graph_as_mol(molecule)
    # molecule = halogenation(molecule, "cl")
    # drawer.draw_graph_as_mol(molecule)
    molecule = halogenation(molecule, "br")
    drawer.draw_graph_as_mol(molecule)


if __name__ == '__main__':
    main()
