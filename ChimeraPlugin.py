import chimera
from chimera import runCommand
from chimera import openModels, Molecule, Vector, selection
import _surface
from numpy import array, single as floatc, intc

class ChimeraPlugin():
    @staticmethod
    def return_object_list():
        mol_list = []
        for mol in openModels.list(modelTypes=[Molecule]):
            mol_list.append(mol.name)
        return mol_list

    def select(self, selection=None):
        if selection is None:
            return runCommand("select")
        else:
            runCommand(str("select ") + selection)

    def save_selection(self, name):
        runCommand(str("namesel " + str(name)))

    def current_atoms(self):
        return selection.currentAtoms()

    @staticmethod
    def save(model_number, path):
        print path
        runCommand("write format pdb #" + model_number + " " + path)

    @staticmethod
    def make_icosahedron(model_number, x, y, z):
        runCommand("shape icosahedron radius 0.5 center " + str(x) + "," + str(y) + "," + str(z) + " modelId #" +
                   str(model_number))

    @staticmethod
    def run(script):
        runCommand("open " + str(script))

    @staticmethod
    def do(what, param = None):
        runCommand(what + str(param))

    @staticmethod
    def append_node(molecule, list, residue):
        at = molecule.newAtom("Tunnel", chimera.Element("Tunn"))
        at.setCoord(chimera.Coord(list[0], list[1], list[2]))
        at.radius = list[3]
        residue.addAtom(at)


    def parse_PDB_channel(self, file_path, name):
        file = open(file_path)
        tunnelObject = chimera.Molecule()
        tunnel = tunnelObject.newResidue(name, " ", 1, " ")
        tunnel_number = 0

        for line in file:

            if line.startswith('HETATM'):

                tunnel_number = int(line[25:30])
                nodes = [float(line[31:39]), float(line[39:47]), float(line[47:55]), float(line[62:68])]
                self.append_node(tunnelObject, nodes, tunnel)

        file.close()
        self.make_channel(tunnelObject, tunnel_number)

    @staticmethod
    def make_channel(model, name):
        model.name = "Tunnel_" + str(name)
        chimera.openModels.add([model])
        chimera.runCommand('repr cpk')

    def get_selection(self, name = None):
        if name is not None:
            chimera.runCommand("select " + name)

        chimera.runCommand("writesel help_file.txt")
        chimera.runCommand("writesel help_model.txt itemType model")

        with open('help_model.txt', 'r') as help_model:

            for mol in openModels.list(modelTypes=[Molecule]):

                for line in help_model:

                    if str(mol) == str(line[:-1]):

                        return mol.name


