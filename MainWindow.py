# coding=utf-8

import os
import sys
from xml.dom.minidom import parse
from importlib.util import spec_from_file_location, module_from_spec
from datetime import datetime
import pickle
import urllib
import time
import random
import socket
import json
import re
import webbrowser

import PointDialog

import Manager

try:
    from PymolPlugin import PymolPlugin as plugin
    PLUGIN = "PyMOL"
    from pymol.Qt import QtCore, QtWidgets
    Qt = QtCore.Qt

except:
    import chimera
    from ChimeraPlugin import ChimeraPlugin as plugin
    PLUGIN = "Chimera"
    from PyQt5 import QtCore, QtWidgets

# Globals:
CONFIGFILE = '.MOLE_PluginSettings.txt'


class MainDialog(QtWidgets.QDialog):
    name = "MOLE 2.5"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.plugin_type = PLUGIN

        self.setWindowTitle('MOLE 2.5')
        #self.setFixedSize(self.size()) # ??
        layout = QtWidgets.QVBoxLayout(self)

        # region Create Frame and NoteBook
        #self.mainframe = Frame(self.parent, width=463, height=623)
        self.mainframe = QtWidgets.QFrame()
        layout.addWidget(self.mainframe)
        #self.mainframe.pack(fill='both', expand=1)

        # TODO: Set wrong executable event
        #self.mainframe.bind('<<WrongExecutable>>',
        #                    lambda e: self.when_error(e, 'Your MOLE 2.0 executable was not found!'))
        #root.bind('<F5>', (lambda event: self.set_structures(
        #    self.input_structure_box)))
        # TODO: Set tooltip balloon = Pmw.Balloon(self.mainframe)

        self.points = {}
        # Binary file, Working directory, csa file
        self.main_parameters = ['', '', '']

        #self.notebook = Pmw.NoteBook(self.mainframe)
        self.notebook = QtWidgets.QTabWidget()
        layout.addWidget(self.notebook)
        #self.notebook.pack(fill='both', expand=1, padx=10, pady=10)
        # endregion

        # region self.mainPage / settings
        self.mainpage = QtWidgets.QWidget()
        mainpage_layout = QtWidgets.QVBoxLayout(self.mainpage)
        self.notebook.addTab(self.mainpage, 'Compute Tunnels')
        self.mainpage.setFocus()

        input_structure_group = QtWidgets.QGroupBox('Specify Input Structure')
        input_structure_group_layout = QtWidgets.QVBoxLayout(input_structure_group)
        mainpage_layout.addWidget(input_structure_group)
        #input_structure_group.pack(fill='both')

        starting_point_group = QtWidgets.QGroupBox('Specify Starting Point')
        mainpage_layout.addWidget(starting_point_group)
        #starting_point_group.pack(fill='both')

        if self.plugin_type == "PyMOL":
            initialized_structs = ('all',) + plugin.return_tuple_objects()
        else:
            initialized_structs = ()

        #self.input_structure_box = Pmw.ScrolledListBox(input_structure_group.interior(), items=initialized_structs,
        #                                               labelpos='nw', listbox_height=4, listbox_selectmode=EXTENDED,)
        self.input_structure_box = QtWidgets.QListWidget()
        self.input_structure_box.addItems(initialized_structs)
        self.input_structure_box.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        input_structure_group_layout.addWidget(self.input_structure_box)

        self.show()
        return
        self.input_structure_box.component('listbox').configure(
            exportselection=0, background='white')
        self.input_structure_box.pack(fill='both', expand=0, padx=10, pady=5)
        balloon.bind(self.input_structure_box,
                     'Select one or more structures in which you want to find channels.')

        if self.plugin_type == "Chimera":
            self.input_structure_box.insert('end', 'all')

            if plugin.return_object_list() is not None:

                for item in plugin.return_object_list():
                    self.input_structure_box.insert('end', item)

        query = Pmw.Group(input_structure_group.interior(),
                          tag_text="Not Active Residues")
        query.pack(fill='both')
        self.query_entry = Entry(query.interior(), width=65)
        self.query_entry.config(background='white')
        self.query_entry.grid(column=0, row=0, padx=10,
                              pady=5, sticky=W + E + N + S)
        self.query_entry.columnconfigure(0, minsize=65)
        self.query_entry.bind(
            "<KeyRelease>", lambda event: root.after(2500, self.validate_query))
        self.query_entry.bind('<Control-KeyRelease-a>',
                              lambda event: self.select_all_query_entry())
        self.query_entry.focus_set()

        balloon.bind(self.query_entry, 'Input not active residues here')

        self.query_label = Label(query.interior(),
                                 text='Select atoms/residues not to be included in the calculation using PatternQuery '
                                      'syntax.', width=68)
        self.query_label.grid(column=0, row=1, padx=10,
                              pady=5, sticky=W + E + N + S, columnspan=2)
        self.query_help = Label(
            query.interior(), text='?', fg="blue", cursor="hand2")
        self.query_help.grid(column=1, row=0, padx=[0, 10], pady=5)
        self.query_help.bind("<Button-1>", lambda event: webbrowser.open_new("https://webchem.ncbr.muni.cz/Wiki"
                                                                             "/PatternQuery:UserManual"))
        balloon.bind(self.query_help, 'PatternQuery Wiki pages.')
        self.is_valid = False

        self.start_points_box = Pmw.ScrolledListBox(starting_point_group.interior(), items=(), labelpos='nw',
                                                    listbox_height=4, listbox_selectmode=EXTENDED)
        self.start_points_box.component('listbox').configure(
            exportselection=0, background='white')
        self.start_points_box.pack(fill='both', expand=0, padx=10, pady=5)
        balloon.bind(self.start_points_box,
                     'Starting point list. If no starting point is specified, \nMOLE plugin will try to find suitable '
                     'starting points automatically.\nOtherwise all selected points will be used.')

        self.button_box1 = Pmw.ButtonBox(starting_point_group.interior())
        self.button_box1.pack(fill='both', expand=0, padx=10, pady=5)
        self.button_box1.add('AddStartingPoint', text='Add Starting Point',
                             command=lambda: self.add_point(self.start_points_box))
        self.button_box1.add('RemoveStartingPoint', text='Remove Starting Point',
                             command=lambda: self.remove_point(self.start_points_box))
        self.button_box1.add('RefreshStructures', text='Refresh Structures',
                             command=lambda: self.set_structures(self.input_structure_box))
        self.button_box1.alignbuttons()

        property_group = Pmw.Group(self.mainpage, tag_pyclass=None)
        property_group.pack(fill='both')

        self.overwrite_results = BooleanVar()
        self.overwrite_results.set(True)
        self.overwrite_results_button = Checkbutton(property_group.interior(), text="Overwrite results",
                                                    variable=self.overwrite_results,
                                                    onvalue=True, offvalue=False)
        self.overwrite_results_button.grid(
            column=0, row=0, padx=10, pady=5, sticky=W + E + N + S)
        balloon.bind(self.overwrite_results_button,
                     'If checked MOLE will overwrite old files in output folder. Otherwise, new folder will be '
                     'created in output folder.')

        self.remove_hydrogens = BooleanVar()
        self.remove_hydrogens_button = Checkbutton(property_group.interior(), text="Remove hydrogens",
                                                   variable=self.remove_hydrogens,
                                                   onvalue=True, offvalue=False)
        self.remove_hydrogens_button.grid(
            column=1, row=0, padx=10, pady=5, sticky=W + E + N + S)
        balloon.bind(self.remove_hydrogens_button,
                     'If checked MOLE will remove all hydrogens from the structure prior to the calculation.')

        self.ignore_het = BooleanVar()
        self.ignore_het_button = Checkbutton(property_group.interior(), text="Ignore HETeroatoms",
                                             variable=self.ignore_het,
                                             onvalue=True, offvalue=False)
        self.ignore_het_button.grid(
            column=2, row=0, padx=10, pady=5, sticky=W + E + N + S)
        balloon.bind(self.ignore_het_button,
                     'If checked MOLE will exclude all HETATM entries prior to the calculation.')

        self.select_working_directory_button = Button(property_group.interior(), text='Save output to:',
                                                      command=self.select_working_directory)
        self.select_working_directory_button.grid(
            row=1, column=0, sticky=W + E + N + S, padx=10, pady=5)
        self.select_working_directory_button.columnconfigure(0, weight=1)
        balloon.bind(self.select_working_directory_button,
                     'Where do you wish to save output from MOLE 2.0 plugin.')

        self.working_directory = Pmw.EntryField(
            property_group.interior(), labelpos='w')
        self.working_directory.component('entry').configure(background='white')
        self.working_directory.grid(
            row=1, column=1, columnspan=2, sticky=W + E + N + S, padx=10, pady=5)
        balloon.bind(self.working_directory,
                     'Where do you wish to save output from MOLE 2.0 plugin.')

        self.generate_csa_selections_button = Button(property_group.interior(), text='Generate CSA sites:',
                                                     command=self.CSA_button_click)
        self.generate_csa_selections_button.grid(
            column=0, row=2, sticky=W + E + N + S, padx=10, pady=5)
        balloon.bind(self.generate_csa_selections_button,
                     'Specify structure by writing its PDB ID and press Generate button.')
        self.structure_for_CSA = Pmw.EntryField(
            property_group.interior(), labelpos='w')
        self.structure_for_CSA.component('entry').configure(background='white')
        self.structure_for_CSA.component('entry').bind(
            '<Return>', lambda event: self.CSA_button_click())
        self.structure_for_CSA.grid(
            column=1, row=2, columnspan=2, sticky=W + E + N + S, padx=10, pady=5)
        balloon.bind(self.structure_for_CSA,
                     'Specify structure by writing its PDB ID and press Generate button.')

        self.select_CSA_button = Button(property_group.interior(
        ), text='Select CSA file:', command=self.select_CSA)
        self.select_CSA_button.grid(
            column=0, row=3, sticky=W + E + N + S, padx=10, pady=5)
        balloon.bind(self.select_CSA_button,
                     'Insert location of CSA.dat file containing CSA database for suggesting active sites as a '
                     'starting points.')

        self.CSA = Pmw.EntryField(property_group.interior(), labelpos='w')
        self.CSA.component('entry').configure(background='white')
        self.CSA.grid(column=1, columnspan=2, row=3,
                      sticky=W + E + N + S, padx=10, pady=5)
        balloon.bind(self.CSA,
                     'Insert location of CSA.dat file containing CSA database for suggesting active sites as a '
                     'starting points.')

        self.select_executable_button = Button(property_group.interior(), text='MOLE 2.5 location:',
                                               command=self.select_executable)
        self.select_executable_button.grid(
            column=0, row=4, sticky=W + E + N + S, padx=10, pady=5)
        balloon.bind(self.select_executable_button,
                     'Select proper path to the MOLE 2.0 command line location')

        self.executable = Pmw.EntryField(
            property_group.interior(), labelpos='w')
        self.executable.component('entry').configure(background='white')
        self.executable.grid(column=1, columnspan=2, row=4,
                             sticky=W + E + N + S, padx=10, pady=5)
        balloon.bind(
            self.executable, 'Select proper path to the MOLE 2.5 command line location.')

        self.compute_tunnels_button = Button(self.mainpage, text='Compute Tunnels', font=("Helvetica", 12, "bold"),
                                             command=lambda b='tunnels': Manager.Manager().construct_params_and_run(self,
                                                                                                                    b))
        self.compute_tunnels_button.pack(
            side=BOTTOM, fill='both', expand=0, padx=10, pady=5)

        # endregion

        # region self.paramsPage / params
        self.params_page = self.notebook.add('Settings')
        cavity_parameters_group = Pmw.Group(
            self.params_page, tag_text='Cavity Parameters')
        cavity_parameters_group.pack(fill='both', expand=4)

        self.probe_radius = Pmw.Counter(cavity_parameters_group.interior(), labelpos='w', label_text='Probe Radius',
                                        entryfield_value='3', increment=0.1,
                                        entryfield_validate={'validator': 'real', 'separator': '.', 'min': 1.4,
                                                             'max': 45},
                                        datatype={'counter': 'real', 'separator': '.'})
        self.probe_radius.component('entry').configure(background='white')
        balloon.bind(self.probe_radius,
                     'Radius used for construction of molecular surface.')

        self.interior_threshold = Pmw.Counter(cavity_parameters_group.interior(), labelpos='w',
                                              label_text='Interior Threshold',
                                              entryfield_value='1.25', increment=0.1,
                                              entryfield_validate={'validator': 'real', 'separator': '.', 'min': 0.8,
                                                                   'max': 45.0},
                                              datatype={'counter': 'real', 'separator': '.'})
        self.interior_threshold.component(
            'entry').configure(background='white')
        balloon.bind(self.interior_threshold,
                     'Lower bound of the tunnel radius.')

        self.surface_cover_radius = Pmw.Counter(cavity_parameters_group.interior(), labelpos='w',
                                                label_text='Surface Cover Radius',
                                                entryfield_value='10', increment=0.1,
                                                entryfield_validate={'validator': 'real', 'separator': '.', 'min': 5.0,
                                                                     'max': 25.0},
                                                datatype={'counter': 'real', 'separator': '.'})
        self.surface_cover_radius.component(
            'entry').configure(background='white')
        balloon.bind(self.surface_cover_radius,
                     'Determines the density of tunnel exits on the molecular surface.')

        self.origin_radius = Pmw.Counter(cavity_parameters_group.interior(), labelpos='w', label_text='Origin Radius',
                                         entryfield_value='5', increment=0.1,
                                         entryfield_validate={'validator': 'real', 'separator': '.', 'min': 0.1,
                                                              'max': 10.0},
                                         datatype={'counter': 'real', 'separator': '.'})
        self.origin_radius.component('entry').configure(background='white')
        balloon.bind(self.origin_radius,
                     'Better starting points are localized within the defined radius from the original starting point.')

        filter_settings = Pmw.Group(
            self.params_page, tag_text='Tunnel Parameters')
        filter_settings.pack(fill='both', expand=4)
        self.bottleneck_radius = Pmw.Counter(filter_settings.interior(), labelpos='w', label_text='Bottleneck Radius',
                                             entryfield_value='1.25', increment=0.1,
                                             entryfield_validate={'validator': 'real', 'separator': '.', 'min': 0.1,
                                                                  'max': 6.0},
                                             datatype={'counter': 'real', 'separator': '.'})
        self.bottleneck_radius.component('entry').configure(background='white')
        balloon.bind(self.bottleneck_radius, 'The minimum radius of a tunnel.')

        self.bottleneck_length = Pmw.Counter(filter_settings.interior(), labelpos='w', label_text='Bottleneck Length',
                                             entryfield_value='3', increment=0.1,
                                             entryfield_validate={'validator': 'real', 'separator': '.', 'min': 0.0,
                                                                  'max': 20.0},
                                             datatype={'counter': 'real', 'separator': '.'})
        self.bottleneck_length.component('entry').configure(background='white')
        balloon.bind(self.bottleneck_length,
                     'Length of a possible profile narrower than the Bottleneck Radius')

        self.cutoff_ratio = Pmw.Counter(filter_settings.interior(), labelpos='w', label_text='Cutoff Ratio',
                                        entryfield_value='0.7', increment=0.05,
                                        entryfield_validate={'validator': 'real', 'separator': '.', 'min': 0.0,
                                                             'max': 1.0},
                                        datatype={'counter': 'real', 'separator': '.'})
        self.cutoff_ratio.component('entry').configure(background='white')
        balloon.bind(self.cutoff_ratio,
                     'Determines maximum similarity of tunnels centerline. \nIf two tunnels are more similar than the '
                     'threshold, the longer is discarded.')

        # region load executable location from settingsfile

        temp_path = os.path.normcase(
            str(os.environ['TEMP']) if 'win32' == str.lower(sys.platform) else '/tmp/')

        if os.path.exists(CONFIGFILE):

            try:
                with open(CONFIGFILE, 'r') as f:
                    self.main_parameters = list(
                        map(lambda x: os.path.normcase(x), pickle.load(f)))
                    self.executable.setvalue(self.main_parameters[0])
                    self.working_directory.setvalue(self.main_parameters[1])
                    self.CSA.setvalue(self.main_parameters[2])

            except:
                self.executable.setvalue('')
                self.working_directory.setvalue(temp_path)
                self.CSA.setvalue('')

        else:
            self.working_directory.setvalue(temp_path)
            self.main_parameters[1] = temp_path
            with open(CONFIGFILE, 'w') as f:
                pickle.dump(self.main_parameters, f)
        # endregion

        self.params = (
            self.probe_radius, self.interior_threshold, self.surface_cover_radius, self.origin_radius,
            self.bottleneck_radius,
            self.bottleneck_length, self.cutoff_ratio)
        Pmw.alignlabels(self.params)

        for counter in self.params:
            counter.pack(fill='both', expand=1, padx=10, pady=5)

        weight_function_group = Pmw.Group(
            self.params_page, tag_text='Weight Function')
        weight_function_group.pack(fill='both', expand=0)

        self.weight_function = StringVar()

        self.voronoi_button = Radiobutton(weight_function_group.interior(), text='Voronoi Scale',
                                          variable=self.weight_function, value="VoronoiScale")
        self.voronoi_button.grid(
            column=0, row=0, padx=30, pady=5, sticky=W + E + N + S)
        balloon.bind(self.voronoi_button,
                     'Recommended for identification of tunnels leading to the buried active sites')

        self.length_radius_button = Radiobutton(weight_function_group.interior(), text='Length + Radius',
                                                variable=self.weight_function, value="LengthAndRadius")
        self.length_radius_button.grid(
            column=1, row=0, padx=30, pady=5, sticky=W + E + N + S)
        balloon.bind(self.length_radius_button,
                     'Old and universal MOLE 2.0 weight function')

        self.length_button = Radiobutton(weight_function_group.interior(), text='Length', variable=self.weight_function,
                                         value="Length")
        self.length_button.grid(column=2, row=0, padx=30,
                                pady=5, sticky=W + E + N + S)
        balloon.bind(self.length_button,
                     'Recommended for identification of transmembrane pores')

        self.voronoi_button.invoke()

        # endregion

        # region self.pathPage / settings
        self.path_page = self.notebook.add('Compute Pores')

        input_structure_group = Pmw.Group(
            self.path_page, tag_text='Specify pores starting points')
        input_structure_group.pack(fill='both')

        self.path_starting_points_box = Pmw.ScrolledListBox(input_structure_group.interior(), items=(), labelpos='nw',
                                                            label_text='Starting Points', listbox_height=4,
                                                            listbox_selectmode=EXTENDED)
        self.path_starting_points_box.component('listbox').configure(
            exportselection=0, background='white')
        self.path_starting_points_box.pack(
            fill='both', expand=0, padx=10, pady=5)
        self.path_starting_points_box.bind(
            '<Delete>', lambda event: self.remove_point(self.path_starting_points_box))
        balloon.bind(self.path_starting_points_box,
                     'Starting point list. Every starting point must have coresponding end point.')

        self.button_box2 = Pmw.ButtonBox(input_structure_group.interior())
        self.button_box2.pack(fill='both', expand=0, padx=10, pady=5)
        self.button_box2.add('AddStartingPoint', text='Add Starting Points',
                             command=lambda: self.add_point(self.path_starting_points_box))
        self.button_box2.add('RemoveStartingPoint', text='Remove Starting Points',
                             command=lambda: self.remove_point(self.path_starting_points_box))
        self.button_box2.alignbuttons()

        input_structure_group = Pmw.Group(
            self.path_page, tag_text='Specify pores end points')
        input_structure_group.pack(fill='both')

        self.path_end_points_box = Pmw.ScrolledListBox(input_structure_group.interior(), items=(), labelpos='nw',
                                                       label_text='End Points', listbox_height=4,
                                                       listbox_selectmode=EXTENDED)
        self.path_end_points_box.component('listbox').configure(
            exportselection=0, background='white')
        self.path_end_points_box.bind(
            '<Delete>', lambda event: self.remove_point(self.path_end_points_box))
        self.path_end_points_box.pack(fill='both', expand=0, padx=10, pady=5)
        balloon.bind(self.path_end_points_box,
                     'End points list. Every end point must have coresponding starting point.')

        self.button_box3 = Pmw.ButtonBox(input_structure_group.interior())
        self.button_box3.pack(fill='both', expand=0, padx=10, pady=5)
        self.button_box3.add('AddEndPoint', text='Add End Points',
                             command=lambda: self.add_point(self.path_end_points_box))
        self.button_box3.add('RemoveEndPoint', text='Remove end points',
                             command=lambda: self.remove_point(self.path_end_points_box))
        self.button_box3.alignbuttons()

        self.compute_pores_button = Button(self.path_page, text='Compute Pores', font=("Helvetica", 12, "bold"),
                                           command=lambda b='pores': Manager.Manager().construct_params_and_run(self, b))
        self.compute_pores_button.pack(fill='both', expand=0, padx=10, pady=5)
        # endregion

        # prev set

        # region READ
        self.read_page = self.notebook.add('Read Channels')
        read_group = Pmw.Group(
            self.read_page, tag_text='Select a file with previously computed MOLE tunnels/pores')
        read_group.pack(fill='both')
        self.open_channels_button = Button(read_group.interior(), text='Open computation results',
                                           command=self.open_channels,
                                           font=("Helvetica", 12, "bold"))
        self.open_channels_button.pack(fill='both', expand=0, padx=10, pady=5)

        # endregion

        # region QuickStartGuide

        self.guide_page = self.notebook.add('Quick Guide')

        guide = Text(self.guide_page, width=72, wrap=WORD)
        guide.config(padx=5, pady=8)
        guide.tag_configure('big', font=('Arial', 12, 'bold'))
        guide.tag_configure('plain_text', font=('Arial', 10))
        guide.tag_configure('smaller', font=('Arial', 9))

        guide.insert(END, 'Plugin description:\n\n', 'big')
        plugin_description = ("The plugin is separated into several tabs The crucial for the calculation are: "
                              "Compute Tunnels, Settings, and Compute Pores. At first, specify location of output "
                              "folder and MOLE 2.5 command line executable. Afterwards, select a structure in the "
                              "'input structure listbox' and starting point from the 'starting point listbox'. "
                              "Additionally, if you provide a path to the CSA database [1], the plugin will suggest "
                              "you the potential starting points.\n\n")
        guide.insert(END, plugin_description, 'plain_text')

        guide.insert(END, "Run:\n\n", 'big')
        run_description = ("After selecting one or more structures and one or more starting points by clicking 'Add Starting Point'. Additional "
                           "search parameters can be adjusted in the Settings tabs. For further info on how to use this "
                           "please refer the included manual or visit our webpages. Whenever you would feel lost just "
                           "hover your cursor above any element in order to get tooltip. For more info and news about the "
                           "MOLE 2.5 visit our webpages.\n")

        guide.insert(END, run_description, 'plain_text')

        improve_text = ("\n"
                        "Also if you would like to make a suggestion on how to "
                        "improve MOLE or send a bug report, please contact authors at webchemistryhelp@gmail.com "
                        "or the author of this extension directly - ondra.balcarek@seznam.cz\n\n")

        guide.insert(END, improve_text, 'plain_text')
        guide.insert(END, "Happy tunneling!    ", 'big')
        team_text = ("Mole development team. http://mole.chemi.muni.cz\n\n"
                     "[1] Furnham,N., Holliday,G.L., De Beer,T.A.P., Jacobsen,J.O.B., Pearson,W.R. and Thornton,J.M. (2014) "
                     "The Catalytic Site Atlas 2.0: Cataloging catalytic sites and residues identified in enzymes. Nucleic "
                     "Acids Res., 42, 1–5.")
        guide.insert(END, team_text, 'smaller')
        guide.config(state=DISABLED)
        

        guide.pack(side=LEFT, fill='both')

        # endregion

        # region AuthorsPage
        self.authors_page = self.notebook.add('Authors')
        authors_text_widget = Text(self.authors_page, width=72, height=17)
        authors_text_widget.config(padx=5, bg='#52A300', pady=8)
        authors_text_widget.tag_configure(
            'plain_text', font=('Arial', 10), foreground='white')
        authors_text_widget.insert(
            END, "If you find this tool useful for your work please cite it as:\n\n", 'plain_text')
        cite_text = ("Sehnal D, Svobodova Varekova R, Berka K, Pravda L, Navratilova V,\nBanas P, Ionescu C-M, Otyepka M, "
                     "Koca J. MOLE2.0: advanced approach for analysis of biomacromolecular channels. Journal of "
                     "Cheminformatics 2013, 5:39., doi:10.1186/1758-2946-5-39.\n\n"
                     "If you were using the web server, which is available at http://mole.upol.cz/ please cite it as:\n\n"
                     "Berka K, Hanak O, Sehnal D, Banas P, Navratilova V, Jaiswal D,Ionescu C-M, Svobodova Varekova R, "
                     "Koca J, Otyepka M.\nMOLEonline 2.0: interactive web-based analysis of biomacromolecular\nchannels. Nucleic "
                     "acids research 2012, 40:W222?7., doi:10.1093/nar/gks363")

        authors_text_widget.insert(END, cite_text, 'plain_text')
        authors_text_widget.config(state=DISABLED)
        authors_text_widget.pack(side=LEFT, fill='both')

        # endregion

        Label(self.mainframe, relief='sunken', anchor=W, justify=LEFT, bg='#52A300', fg='white', font=("Helvetica", 12),
              padx=10, pady=10,
              text="(c) 2017 CEITEC & NCBR MU & FCH UPOL\nhttp://mole.chemi.muni.cz                                v. "
                   "2.5.17.7.11").pack(
            fill='both')

        self.notebook.setnaturalsize()
        self.wd = os.path.realpath(os.path.normcase(self.main_parameters[1]))
        self.xml_wd = os.path.realpath(
            os.path.normcase(self.main_parameters[1] + "/xml/"))
        self.profile_wd = os.path.realpath(
            os.path.normcase(self.main_parameters[1] + "/pymol/"))
        os.system('xset r off')  # for keypress

        self.original_view = None

    def keyPressEvent(self, event: QtCore.QEvent):
        if event.key() == Qt.Key_F5:
            self.set_structures(self.input_structure_box)

    def add_point(self, listbox):
        """
        Run add_point window.
        :param listbox:
        :return:
        """
        dialog = PointDialog.PointDialog(
            self, self.mainframe, CSA=self.CSA.get())
        dialog.activate(geometry='centerscreenalways')

        if dialog.return_value is not None:
            r = dialog.point_name + "  Type: " + \
                dialog.return_value[1]['Element']
            x = listbox.get()
            x += (r,)
            listbox.setlist(x)

            self.points[dialog.point_name] = dialog.return_value

    def select_all_query_entry(self, event=None):
        self.query_entry.select_range(0, 'end')
        self.query_entry.icursor('end')

    def remove_point(self, listbox):
        """
        Remove point from listbox
        :param listbox:
        :return:
        """
        x = ()

        for i in listbox.get():

            if i not in listbox.getvalue():
                x += (i,)

        listbox.setlist(x)

        p = {}

        for i in self.points.keys():

            if i not in listbox.getvalue() or i[4:5] == '|':
                p[i] = self.points[i]

        self.points = p

    def select_working_directory(self):
        """
        Open window for selection wd
        :return:
        """
        file_path = tkFileDialog.askdirectory(title='Select Output Directory')
        file_path = os.path.normcase(file_path)

        if file_path != "" and os.access(file_path, os.R_OK and os.W_OK):

            try:
                with open(os.path.join(file_path, 'testfile.txt'), 'w') as f:
                    pass
                with open(os.path.join(file_path, 'testfile.txt'), 'r') as f:
                    pass

            except:
                tkMessageBox.showinfo('Info', 'Selected output directory does not have sufficient permissions to be '
                                              'used.')
                return

            os.remove(os.path.join(file_path, 'testfile.txt'))
            self.working_directory.setvalue(file_path)

            try:
                self.main_parameters[1] = file_path
                with open(CONFIGFILE, 'w') as f:
                    pickle.dump(self.main_parameters, f)

            except:
                pass
        else:
            tkMessageBox.showinfo(
                'Info', 'Selected output directory does not have sufficient permissions to be used.')

    def select_CSA(self):
        """
        Open dialog for CSA file selection.
        :return:
        """
        file_path = tkFileDialog.askopenfilename(title='Select CSA.dat file',
                                                 filetypes=[('CSA.dat file', '.dat'), ('.txt file', '.txt')])
        file_path = os.path.normcase(file_path)

        if file_path != "":
            self.CSA.setvalue(file_path)

            try:
                self.main_parameters[2] = file_path
                with open(CONFIGFILE, 'w') as f:
                    pickle.dump(self.main_parameters, f)

            except:
                pass

    def select_executable(self):
        """
        Open dialog for executable selection.
        :return:
        """
        file_path = tkFileDialog.askopenfilename(title='Select MOLE 2.0 Executable',
                                                 filetypes=[('MOLE 2.0 Executable', '.exe')])
        file_path = os.path.normcase(file_path)

        if file_path != "":
            self.executable.setvalue(file_path)

            try:
                self.main_parameters[0] = file_path
                with open(CONFIGFILE, 'w') as f:
                    pickle.dump(self.main_parameters, f)

            except:
                pass

    def rename_buttons(self, param):
        """
        Rename buttons from tunnels to pores, if necessary.
        :param param:
        :return:
        """
        if param == 'tunnels':
            self.compute_tunnels_button.config(text='Compute Tunnels')

        else:
            self.compute_pores_button.config(text='Compute Pores')

    def CSA_button_click(self):
        """
        Calls CSA selection.
        :return:
        """
        Manager.Manager().generate_CSA_selections(self)

    def set_state(self, state):
        """
        Change state of elements in GUI.
        :param state:
        :return:
        """
        self.generate_csa_selections_button.config(state=state)
        self.structure_for_CSA.component('entry').config(state=state)
        self.ignore_het_button.config(state=state)
        self.overwrite_results_button.config(state=state)
        self.compute_tunnels_button.config(state=state)
        self.compute_pores_button.config(state=state)
        self.select_executable_button.config(state=state)
        self.select_working_directory_button.config(state=state)
        self.select_CSA_button.config(state=state)
        self.CSA.component('entry').config(state=state)
        self.executable.component('entry').config(state=state)
        self.working_directory.component('entry').config(state=state)
        self.remove_hydrogens_button.config(state=state)

        for w in self.params:
            w.component('entry').config(state=state)

        self.button_box1.button(0).config(state=state)
        self.button_box1.button(1).config(state=state)
        self.button_box1.button(2).config(state=state)
        self.button_box2.button(0).config(state=state)
        self.button_box2.button(1).config(state=state)
        self.button_box3.button(0).config(state=state)
        self.button_box3.button(1).config(state=state)

    def when_error(self, event, message):
        """
        Change state to normal and open tkmessagebox.
        :param event:
        :param message:
        :return:
        """
        self.set_state('normal')
        tkMessageBox.showinfo('Error', message)

    def set_structures(self, listbox):
        """
        Fill input structures listbox
        :param listbox:
        :return:
        """
        listbox.clear()
        listbox.insert('end', 'all')

        for item in plugin.return_object_list():
            listbox.insert('end', item)

    def get_element(self, selection):
        """
        Get element from point structure
        :param selection:
        :return:
        """
        for key in self.points.keys():
            if key in selection:

                return self.points[key]

        return '42'

    def open_channels(self):
        file_path = tkFileDialog.askopenfilenames(title='Select MOLE tunnel files',
                                                  filetypes=[('MOLE tunnels', '.xml .pdb .py')])

        if self.plugin_type == "PyMOL":
            file_path = self.fix_list(file_path)

        for chan_file in file_path:

            chan_file = os.path.normcase(chan_file)
            extension = os.path.splitext(chan_file)[1].lower()
            name = os.path.basename(chan_file).split('.')[0].lower()

            if len(extension) < 3:
                return

            try:
                if extension == '.py':
                    plugin.run(chan_file)
                if extension == '.xml':
                    Manager.Manager().parse_XML_channel(chan_file, name, self)
                if extension == '.pdb':
                    plugin().parse_PDB_channel(chan_file, name)

            except Exception as e:
                print(e)

    @staticmethod
    def fix_list(filenames):
        # do nothing if already a python list
        if isinstance(filenames, list):
            return filenames

        # http://docs.python.org/library/re.html
        # the re should match: {text and white space in brackets} AND anynonwhitespacetokens
        # *? is a non-greedy match for any character sequence
        # \S is non white space

        # split filenames string up into a proper python list
        result = re.findall(r"{.*?}|\S+", filenames)

        # remove any {} characters from the start and end of the file names
        result = [re.sub("^{|}$", "", i) for i in result]

        return result

    def group_objects(self, script_path, script_name, structure):
        """

        :param script_path:
        :param script_name:
        :param structure:
        :return:
        """
        if os.path.exists(script_path):

            previous = plugin.return_object_list()
            spec = spec_from_file_location(
                structure + str(datetime.now().toordinal()), script_path)
            module = module_from_spec(spec)
            spec.loader.exec_module(module)
            #imp.load_source(
            #    structure + str(datetime.now().toordinal()), script_path)

            for o in plugin.return_object_list():

                if o not in previous:
                    plugin.return_group(script_name, o)

    def when_computation_done(self, out, param):
        """
        After computation of tunnels is done, run tunnel script
        :param out:
        :return:
        """
        tkMessageBox.showinfo('Computation is done', out)
        self.set_state('normal')

        if self.plugin_type == "Chimera":
            self.profile_wd = os.path.realpath(
                os.path.normcase(self.main_parameters[1] + "/chimera/"))

        tunnels_script = os.path.join(self.profile_wd, 'complex.py')
        paths_script = os.path.join(self.wd, 'paths.py')
        pores_script = os.path.join(self.wd, 'pores.py')
        autopores_script = os.path.join(self.wd, 'autopores.py')
        dom = parse(os.path.join(self.xml_wd, 'tunnels.xml'))

        for e in dom.getElementsByTagName('Exception'):
            tkMessageBox.showinfo('Exception', e.getAttribute('Text'))

        if self.plugin_type == "PyMOL":

            if param == "pores":
                self.group_objects(tunnels_script, 'Pores', "MOLE_tunnels")
            else:
                self.group_objects(tunnels_script, 'Tunnels', "MOLE_pores")

            self.group_objects(paths_script, 'Paths', 'MOLE_Paths')
            self.group_objects(pores_script, 'Pores', 'MOLE_Pores')
            self.group_objects(autopores_script, 'AutoPores', 'MOLE_AutoPores')
            plugin().show("cgo", "Interior*")
            plugin().set_view(self.original_view)

        else:
            plugin.run(tunnels_script)

        self.set_state('normal')

    def validate_query(self, event=None):
        """
        Checks internet connection
        :param pdb_id:
        :return: validation report
        """

        if len(self.query_entry.get()) == 0:
            return

        query = self.strip_white(self.query_entry.get())
        query_page = 'http://webchem.ncbr.muni.cz/Platform/PatternQuery/ValidateQuery?query=' + query
        query_page = urllib.parse.quote(query_page, ':/' + '?' + '=')
        data_response = False
        tries = 0
        limit = 5

        while tries < limit:

            try:
                req = urllib.request.Request(query_page)
                response = urllib.request.urlopen(req, None, 5)

            except urllib.error.HTTPError as e:
                if e.code == 404:
                    self.when_error('', "Could not find server!")

                else:
                    time.sleep(random.randint(1, 2))  # sleep 1-2 seconds
                    tries += 1

            except urllib.error.URLError as e:
                time.sleep(random.randint(1, 2))  # sleep 1-2 seconds
                tries += 1

            except socket.timeout as e:
                time.sleep(random.randint(1, 2))  # sleep 1-2 seconds
                tries += 1

            else:
                data_response = True
                self.isOk = json.load(response)
                break

        if data_response == False:
            self.when_error('', "Check your internet connection!")
            return

        if 'error' in self.isOk:
            self.query_label.config(
                text="Invalid query! Error: " + self.isOk['error'])
            self.is_valid = False
            return

        else:
            self.query_label.config(text="Valid!")

            self.is_valid = True
            return

    def fill_boxes(self, content):
        content1 = tuple(
            filter(lambda x: x not in self.start_points_box.get(), content))
        content2 = tuple(
            filter(lambda x: x not in self.path_starting_points_box.get(), content))
        content3 = tuple(
            filter(lambda x: x not in self.path_end_points_box.get(), content))

        self.start_points_box.setlist(
            tuple(self.start_points_box.get()) + content1)
        self.path_starting_points_box.setlist(
            tuple(self.path_starting_points_box.get()) + content2)
        self.path_end_points_box.setlist(
            tuple(self.path_end_points_box.get()) + content3)

    def strip_white(self, text):
        lst = text.split('"')
        for i, item in enumerate(lst):
            if not i % 2:
                lst[i] = re.sub(r"\s+", "", item)
        return str('"'.join(lst))

    def __del__(self):
        os.system('xset r on')


if PLUGIN == "Chimera":
    chimera.dialogs.register(MainWindow.name, MainWindow)
    dir, file = os.path.split(__file__)
    icon = os.path.join(dir, 'molelogo.png')
    chimera.tkgui.app.toolbar.add(
        icon, lambda d=chimera.dialogs.display, n=MainWindow.name: d(n), 'MOLE 2.5', None)

########## DEBUG ##########

if __name__ == '__main__':
    app = QtWidgets.QApplication([])
    window = QtWidgets.QMainWindow()
    app.setActiveWindow(window)
    widget = MainDialog(window)
    app.exec_()
