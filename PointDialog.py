import Pmw
from Tkinter import *
import tkMessageBox
try:
    import PymolPlugin
    PLUGIN = "PyMOL"
    from pymol.cgo import *

except(NameError, ImportError):
    import ChimeraPlugin
    PLUGIN = "Chimera"
    import re


class PointDialog(Pmw.Dialog):
    """Creates point dialog for adding starting or end point"""

    def __init__(self, mole_object,parent=None, **kw):

        # Initialise base class (after defining options).
        Pmw.Dialog.__init__(self, parent, command=self.return_to_mainwindow)

        self.plugin_type = PLUGIN
        self.CSA = ''
        if kw.has_key('CSA'):
            self.CSA = kw['CSA']

        self.point_name = 'Point007'

        if kw.has_key('PointName'):
            self.point_name = kw['PointName']

        self.structure = ''
        if kw.has_key('Structure'):
            self.structure = kw['Structure']

        # Create the components.
        interior = self.interior()

        balloon = Pmw.Balloon(interior)

        self.point_group = []
        self.point_var = IntVar()
        self.point_var.set(0)

        # region point selection
        point_box = Pmw.Group(interior, tag_pyclass=Radiobutton, tag_text=mole_object.plugin_type + ' selection (eg.: \"(sele)\")',
                      tag_value=0, tag_variable=self.point_var)
        point_box.pack(fill='both', expand=1)

        self.point_0 = Pmw.EntryField(point_box.interior(), labelpos='w', label_text='Point:', value='(sele)')
        self.point_0.pack(fill='both', expand=1, padx=10, pady=5)
        balloon.bind(self.point_0, 'Example: (sele)\nThe point is defined as a center of mass of selected residues. \n'
                     + 'You can choose residues using by selecting them in PyMOL and writing selection name. \n')

        self.point_group.append(point_box)
        # endregion

        # region starting point coords
        point_box = Pmw.Group(interior, tag_pyclass=Radiobutton, tag_text='Point coordinates',
                      tag_value=2, tag_variable=self.point_var)
        point_box.pack(fill='both', expand=1)

        self.point_x = Pmw.Counter(point_box.interior(), labelpos='w', label_text='X',
                                   entryfield_value=0, increment=0.1,
                                   entryfield_validate={'validator': 'real', 'separator': '.'},
                                   datatype={'counter': 'real', 'separator': '.'},
                                   entryfield_command=self.show_crisscross)
        self.point_x.component('uparrow').bind("<Button 1>", self.when_arrow_clicked, add="+")
        self.point_x.component('downarrow').bind("<Button 1>", self.when_arrow_clicked, add="+")

        self.point_y = Pmw.Counter(point_box.interior(), labelpos='w', label_text='Y',
                                   entryfield_value=0, increment=0.1,
                                   entryfield_validate={'validator': 'real', 'separator': '.'},
                                   datatype={'counter': 'real', 'separator': '.'},
                                   entryfield_command=self.show_crisscross)
        self.point_y.component('uparrow').bind("<Button 1>", self.when_arrow_clicked, add="+")
        self.point_y.component('downarrow').bind("<Button 1>", self.when_arrow_clicked, add="+")

        self.point_z = Pmw.Counter(point_box.interior(), labelpos='w', label_text='Z',
                                   entryfield_value=0, increment=0.1,
                                   entryfield_validate={'validator': 'real', 'separator': '.'},
                                   datatype={'counter': 'real', 'separator': '.'},
                                   entryfield_command=self.show_crisscross)
        self.point_z.component('uparrow').bind("<Button 1>", self.when_arrow_clicked, add="+")
        self.point_z.component('downarrow').bind("<Button 1>", self.when_arrow_clicked, add="+")

        self.compute_center_button = Button(interior, text='Compute Center',
                                            command=lambda :self.compute_center(mole_object))

        c = (self.point_x, self.point_y, self.point_z, self.compute_center_button)
        Pmw.alignlabels(c)

        for i in c:
            i.pack(fill='both', expand=1, padx=10, pady=5)

        self.point_group.append(point_box)
        # endregion

        Pmw.aligngrouptags(self.point_group)

        # Check keywords and initialise options.
        self.initialiseoptions()


    def compute_center(self, mole_object):
        """
        Compute center from selection.
        :param mole_object:
        :return:
        """
        if mole_object.plugin_type == "PyMOL":
            sel = PymolPlugin.PymolPlugin().get_model('all')
            cnt = len(sel.atom)

        else:
            sel = ChimeraPlugin.ChimeraPlugin().select()
            cnt = len(ChimeraPlugin.ChimeraPlugin().current_atoms())

        cent_x = 0
        cent_y = 0
        cent_z = 0

        if cnt == 0:
            return 0, 0, 0

        if mole_object.plugin_type == "PyMOL":

            for a in sel.atom:
                cent_x += a.coord[0]
                cent_y += a.coord[1]
                cent_z += a.coord[2]

        else:

            for a in ChimeraPlugin.ChimeraPlugin().current_atoms():
                cent_x += a.coord()[0]
                cent_y += a.coord()[1]
                cent_z += a.coord()[2]

        cent_x /= cnt
        cent_y /= cnt
        cent_z /= cnt

        self.point_x.component('entryfield').setentry(cent_x)
        self.point_y.component('entryfield').setentry(cent_y)
        self.point_z.component('entryfield').setentry(cent_z)

        self.show_crisscross(mole_object)

    def show_crisscross(self, mole_object):
        """
        Show center of comutation.
        :return:
        """
        if mole_object.plugin_type == "PyMOL":
            obj = [
                LINEWIDTH, 3,

                BEGIN, LINE_STRIP,
                VERTEX, float(float(self.point_x.get()) - 0.5), float(self.point_y.get()), float(self.point_z.get()),
                VERTEX, float(float(self.point_x.get()) + 0.5), float(self.point_y.get()), float(self.point_z.get()),
                END,

                BEGIN, LINE_STRIP,
                VERTEX, float(self.point_x.get()), float(float(self.point_y.get()) - 0.5), float(self.point_z.get()),
                VERTEX, float(self.point_x.get()), float(float(self.point_y.get()) + 0.5), float(self.point_z.get()),
                END,

                BEGIN, LINE_STRIP,
                VERTEX, float(self.point_x.get()), float(self.point_y.get()), float(float(self.point_z.get()) - 0.5),
                VERTEX, float(self.point_x.get()), float(self.point_y.get()), float(float(self.point_z.get()) + 0.5),
                END

            ]

            PymolPlugin.PymolPlugin().delete(self.point_name)
            view = PymolPlugin.PymolPlugin().get_view()
            PymolPlugin.PymolPlugin().load_CGO(obj, self.point_name)
            PymolPlugin.PymolPlugin().set_view(view)

        else:
            chimera_model_number = int(mole_object.input_structure_box.index('active')) - 1
            ChimeraPlugin.ChimeraPlugin().make_icosahedron(str(chimera_model_number), float(self.point_x.get()),
                                                           float(self.point_y.get()), float(self.point_z.get()))

    def when_arrow_clicked(self, event):
        self.show_crisscross()

    def return_to_mainwindow(self, result):
        self.return_value = None

        if result == 'OK' and self.point_var.get() == 0:


            if self.point_0.get() == '':

                self.when_error('', 'Such selection does not exist in PyMOL.')
                return

            try:
                if self.plugin_type == "PyMOL":
                    PymolPlugin.PymolPlugin().get_model(self.point_0.get())

                else:
                    ChimeraPlugin.ChimeraPlugin().get_selection(self.point_0.get())



            except:
                self.when_error('', 'Point Selection is not defined!')
                return

            self.return_value = self.compute_return_value(self.point_var.get())

            if self.return_value[1]['Element'] == 'Residue':
                self.point_name = self.point_0.get()

            else:
                self.point_name = '[' + str(round(float(self.return_value[1]['X']), 3)) + ', ' + str(
                    round(float(self.return_value[1]['Y']), 3)) + ', ' + str(
                    round(float(self.return_value[1]['Z']), 3)) + ']'
        self.deactivate()

    @staticmethod
    def when_error(event, message):
        tkMessageBox.showinfo('Error', message)

    def compute_return_value(self, type):
        value = ()

        if type == 0:


            #try:
            if self.plugin_type == "PyMOL":
                value += ({'Structure': self.structure, 'Type': 'Selection_' + self.point_0.get()},)
                model = PymolPlugin.PymolPlugin().get_model(self.point_0.get())
                help = []

                for atom in model.atom:
                    help += [[atom.chain, atom.resn, str(atom.resi)]]

                origin = self.distinct(help)
                for a in origin:
                    value += ({'Element': 'Residue',
                               'Name': a[1],
                               'Chain': a[0],
                               'SequenceNumber': a[2]},)

                return value

            else:
                value += ({'Structure': ChimeraPlugin.ChimeraPlugin().get_selection(),
                           'Type': 'Selection_' + self.point_0.get()},)
                sequence_number = []

                with open('help_file.txt', 'r') as help_file:

                    for line in help_file:

                        if line.startswith("#"):
                            line = line[3:]

                        start = line[0:3]
                        end = line[-2:]
                        help_sequence = re.search('%s(.*)%s' % (start, end), line).group(1)
                        help_sequence = help_sequence[1:-1]
                        sequence_number.append(help_sequence)

                        value += ({'Element': 'Residue',
                                   'Name': line[0:3],
                                   'Chain': line[-2:-1],
                                   'SequenceNumber': sequence_number[0]}, )

                return value


            #except:
                #self.when_error('', 'Starting Point Selection is not defined!')
                #return

        if type == 2:
            value += ({'Structure': self.structure, 'Type': 'Point'},)
            value += ({'Element': 'Point',
                       'X': self.point_x.get(),
                       'Y': self.point_y.get(),
                       'Z': self.point_z.get()},)

        return value

    @staticmethod
    def distinct(l):
        distinct = []

        for i in l:

            if i not in distinct:
                distinct += [i]

        return distinct
