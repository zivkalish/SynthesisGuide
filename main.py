import drawer, builder
from rdkit.Chem import MolFromSmiles

def main():
    s = "CCC"
    mol = MolFromSmiles(s)
    drawer.draw_mol(mol)
    # graph = builder.molecule2graph(mol)
    # drawer.draw_graph(graph)

if __name__ == '__main__':
    main()
