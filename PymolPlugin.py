from pymol import cmd, selector
from pymol.cmd import _feedback, fb_module, fb_mask, is_list, _cmd
from pymol.cgo import *
from chempy import Bond, Atom
from chempy.models import Indexed
from chempy import Bond, Atom




class PymolPlugin:
    @staticmethod
    def return_tuple_objects():
        """

        :return: Object list as tuple from pymol
        """
        return tuple(cmd.get_object_list())

    @staticmethod
    def return_object_list(s=None):
        if s is None:
            return cmd.get_object_list()
        if s is not None:
            return cmd.get_object_list(s)

    @staticmethod
    def select_structure(csa_struct, y, x):
        """

        :param csa_struct:
        :param y:
        :param x:
        :return: selection of defined params
        """
        return cmd.select(csa_struct + '_' + y, x)

    @staticmethod
    def return_group(name, count):
        """

        :param name: name of group
        :param count: counting
        :return: group
        """
        return cmd.group(name, count)

    @staticmethod
    def get_view():
        """

        :return:
        """
        return cmd.get_view()

    @staticmethod
    def do(what, param):
        return cmd.do(what + param)

    @staticmethod
    def center():
        cmd.center()

    @staticmethod
    def set_view(view):
        """

        :return:
        """
        return cmd.set_view(view)

    @staticmethod
    def get_model(model):
        return cmd.get_model(model)

    @staticmethod
    def delete(what):
        return cmd.delete(what)

    @staticmethod
    def run(what):
        return cmd.run(what)

    @staticmethod
    def save(path, param):
        return cmd.save(path, '"' + param + '"')

    @staticmethod
    def load_CGO(object, point):
        return cmd.load_cgo(object, point)

    @staticmethod
    def show(representation, model):
        return cmd.show(representation, model)

    @staticmethod
    def make_channel(model, i, name=None):

        for a in range(len(model.atom) - 1):
            bd = Bond()
            bd.index = [a, a + 1]
            model.bond.append(bd)

        if name is None:
            name = "Tunnel" + str(i)

        cmd.load_model(model, name, state=1)
        cmd.set("sphere_mode", "0", name)
        cmd.set("sphere_color", "red", name)

        cmd.show("spheres", name)
        cmd.group("Tunnels", name)

        return Indexed()

    @staticmethod
    def append_node(model, list, test = None):
        at = Atom()
        at.name = '0'
        at.vdw = list[3]
        at.coord = list[:3]
        model.atom.append(at)


    def parse_PDB_channel(self, file_path, name):
        file = open(file_path)
        i = 0
        channelId = 1
        model = Indexed()
        nodes = []

        for line in file:

            if line.startswith('HETATM'):

                if channelId != int(line[25:30]):
                    i += 1
                    model = self.make_channel(model, i, name)
                    channelId += 1

                nodes = [float(line[31:39]), float(line[39:47]), float(line[47:55]), float(line[62:68])]
                self.append_node(model, nodes)

        file.close()
        i += 1
        self.make_channel(model, i, name)


