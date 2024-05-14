import os
import subprocess
import sys
import datetime
import threading
import string

from datetime import datetime
from xml.dom.minidom import Document
try:
   from chempy.models import Indexed
except:
    from chimera import Molecule
    pass
from xml.dom.minidom import parse
import csv

try:
    from PymolPlugin import PymolPlugin as plugin

except (ImportError, NameError):
    from ChimeraPlugin import ChimeraPlugin as plugin

# GLOBALS#
OUTPUT = ''


class Manager():
    def run_sub_process(self, path, param):
        """
        Run executable file
        :param path:
        :return:
        """

        if os.path.exists(os.path.realpath(path)):
            plat = str.lower(sys.platform)
            args = [os.path.realpath(path), os.path.join(self.wd, 'params.xml')]

            output = ''
            if 'win32' in plat:
                output = subprocess.Popen(args, stdout=subprocess.PIPE).communicate()[0].split('\n')

            else:
                output = subprocess.Popen(['mono'] + args, stdout=subprocess.PIPE).communicate()[0].split('\n')

            if param == 'tunnels':
                output = [e for e in output if e.startswith("Computing tunnels")][0]
            else:
                output = [e for e in output if e.startswith("Computing pores")][0]
            split_out = output.split()
            output = 'Found ' + split_out[2] + ' ' + param + '.'

            if 0 == int(split_out[2]):
                # output += '\n\nTry adjusting settings in the \'Settings\' tab.'
                output = 'Found ' + split_out[2] + ' ' + param + '\n\nTry adjusting settings in the \'Settings\' tab.'


        else:
            output = 'Your MOLE 2.0 executable was not found!'

        global OUTPUT
        OUTPUT = output


    @staticmethod
    def parse_XML_channel(file_path, name, mole_object):
        xml = parse(file_path)
        pool = ['Tunnel', 'Pore', 'Path']

        for element in pool:
            i = 1

            for tunnel in xml.getElementsByTagName(element):

                if mole_object.plugin_type == "PyMOL":
                    model = Indexed()
                    individual_model = None

                else:
                    model = Molecule()
                    individual_model = model.newResidue("Tunnel_" + str(i), " ", 1, " ")

                for node in tunnel.getElementsByTagName('Node'):

                    nodes = [float(node.getAttribute('X')), float(node.getAttribute('Y')),
                             float(node.getAttribute('Z')), float(node.getAttribute('Radius'))]

                    plugin.append_node(model, nodes, individual_model)

                model = plugin.make_channel(model, i)
                i += 1

        if mole_object.plugin_type == "PyMOL":
            plugin.center()


    def generate_CSA_selections(self, mole_object):
        active_sites = ()
        csa_structure = mole_object.structure_for_CSA.get().lower()

        if mole_object.main_parameters[2] == '':
            mole_object.when_error('', 'Specify CSA.dat file for search of active sites.')
            return

        if csa_structure == '':

            if len(plugin.return_object_list()) == 0:
                mole_object.when_error('', 'No structure was found in object list.')
                return


            mole_object.when_error('', 'No structure was defined.')
            return

        else:
            if csa_structure not in map(lambda x: x.lower(), plugin.return_object_list()):
                mole_object.when_error('', csa_structure + ' was not found in object list.')
                return


        try:
            file = csv.reader(open(mole_object.main_parameters[2], 'r'))
            sites = 0
            s = ''
            listbox_entry = ()
            points_entry = ()
            points_entry += ({'Structure': '', 'Type': 'CSA'},)

            data = []

            for row in file:
                if row[0] == csa_structure:
                    data.append(row)
                    sites = long(row[1])

            for i in range(sites + 1):

                for row in data:

                    if row[1] == str(i):
                        s += ' (chain ' + row[3] + ' & resi ' + row[
                            4] + ') |'  # (row[3] + '/' + row[2] + ' ' + row[4] + '/',)
                        points_entry += ({'Element': 'Residue',
                                          'Name': row[2],
                                          'Chain': row[3],
                                          'SequenceNumber': row[4]},)
                        listbox_entry += (row[3] + row[4],)

                if len(s) > 1:
                    y = ''

                    for j in range(len(listbox_entry) - 1):
                        y += listbox_entry[j] + ','

                    y += listbox_entry[len(listbox_entry) - 1] + ','
                    x = csa_structure + ' & (' + s[0:-1] + ')'

                    if mole_object.plugin_type == "PyMOL":
                        plugin.select_structure(csa_structure, y, x)

                    else:
                        chimera_active_site = y.split(',')
                        del chimera_active_site[-1]
                        swapped = []
                        selection_list = []
                        chimera_sel = "#" + str(mole_object.input_structure_box.index('active') - 1) + ":"

                        for i in chimera_active_site:
                            swapped.append(i[-1] + i[1:-1] + i[0])

                        for i in range(len(swapped)):
                            selection_list.append(swapped[i][:-1] + "." + swapped[i][-1:] + ":")

                        for i in selection_list:
                            chimera_sel += i

                        chimera_sel = chimera_sel[:-1]
                        plugin().select(chimera_sel)
                        plugin().save_selection(x)

                    mole_object.points[csa_structure + '|' + y[0:-1]] = points_entry
                    active_sites += (y,)
                    s = ''
                    listbox_entry = ()
                    points_entry = ()
                    points_entry += ({'Structure': '', 'Type': 'CSA'},)

            mole_object.fill_boxes(tuple((map(lambda x: csa_structure + '|' + x, active_sites))))
            mole_object.structure_for_CSA.clear()

            if len(active_sites) < 1:
                mole_object.when_error('', 'No active sites found for the \'' + csa_structure + '\'')

        except Exception as e:
            mole_object.when_error('',
                                   'An error occurred during the processing of CSA file. Did you provide it in a '
                                   'correct format?')
            print(e)
            return

    def construct_params_and_run(self, mole_object, param):
        """
        Generates params.xml file
        :param mole_object:
        :param param:
        :return:
        """
        mole_object.set_state('disabled')

        if mole_object.plugin_type == "PyMOL":
            mole_object.original_view = plugin.get_view()

        tunnelParameters = {}
        cavityParameters = {}

        for i in mole_object.params:
            tunnelParameters[i['label_text'].replace(' ', '')] = i.get()

        cavityParameters["InteriorThreshold"] = tunnelParameters["InteriorThreshold"]
        cavityParameters["ProbeRadius"] = tunnelParameters["ProbeRadius"]

        del tunnelParameters["ProbeRadius"]
        del tunnelParameters["InteriorThreshold"]

        tunnelParameters['BottleneckTolerance'] = tunnelParameters.pop('BottleneckLength')
        tunnelParameters['MaxTunnelSimilarity'] = tunnelParameters.pop('CutoffRatio')

        if len(plugin.return_object_list()) == 0:
            mole_object.when_error('', 'No structure loaded!')

        if not os.path.exists(mole_object.main_parameters[0]):
            mole_object.when_error('', 'Your MOLE 2.5 executable was not found!')
            return

        # region create working directory and test permissions
        if not os.path.exists(mole_object.main_parameters[1]):
            mole_object.when_error('', 'No output directory specified.')
            return
        # endregion

        self.wd = os.path.realpath(os.path.normcase(mole_object.main_parameters[1]))

        # region create subdirectory
        if not mole_object.overwrite_results.get():
            self.wd = os.path.join(self.wd, str(datetime.datetime.now().toordinal()))

        if not os.path.exists(self.wd):

            try:
                os.mkdir(self.wd)

            except:
                mole_object.when_error('', 'Output directory cannot be created!')
                return

        used_struct = 'All'
        list_of_structures = []

        doc = Document()
        root = doc.createElement("Tunnels")
        doc.appendChild(root)

        if len(mole_object.input_structure_box.getvalue()) < 1:
            mole_object.when_error('', 'No structure selected')
            return

        elif 'all' in mole_object.input_structure_box.getvalue():

            if mole_object.plugin_type == "PyMOL":

                for i in plugin.return_object_list():
                    list_of_structures.append(i)
            else:
                for i in mole_object.input_structure_box.get(1, 'end'):
                    list_of_structures.append(i)

        else:
            selection = set(mole_object.input_structure_box.getvalue())

            for s in selection:
                if mole_object.plugin_type == "PyMOL":
                    current_structure = plugin.return_object_list(s)

                else:
                    chimera_index = int(mole_object.input_structure_box.index('active')) - 1
                    chimera_index = str("#") + str(chimera_index)
                    current_structure = list((mole_object.input_structure_box.getvalue()))

                if current_structure is None:
                    mole_object.when_error('',
                                           'The structure you have selected is no longer available in the PyMOL object '
                                           'list. Please refresh structures.')
                    return

                else:
                    list_of_structures.append(current_structure[0])

            if len(list_of_structures) == 1:
                used_struct = str(list_of_structures[0])

        for i in list_of_structures:
            used_struct = str(i)

            # region Input
            path = os.path.realpath(os.path.join(self.wd, used_struct + '.pdb'))
            e = doc.createElement("Input")

            if len(list_of_structures) == 1:
                if mole_object.plugin_type == "PyMOL":
                    path = os.path.realpath(os.path.join(self.wd, used_struct + '.pdb'))
                    plugin.save(path, list_of_structures[0])

                if mole_object.plugin_type == "Chimera":
                    name = self.wd + "\\" + used_struct + ".pdb"
                    plugin.save(str(mole_object.input_structure_box.index('active') - 1), name)

            else:
                ex = ''
                for i in list_of_structures:
                    ex += i + "|"

                ex = string.rstrip(ex, '|')

                if mole_object.plugin_type == "PyMOL":
                    plugin.do("select x,", ex)
                    plugin.save(path, "x")
                    plugin.delete("x")

            e.appendChild(doc.createTextNode(path))
            root.appendChild(e)
            # endregion

            # region WorkingDirectory
            e = doc.createElement("WorkingDirectory")
            e.appendChild(doc.createTextNode(self.wd))
            root.appendChild(e)
            # endregion

            # region Params
            query_str = mole_object.strip_white(mole_object.query_entry.get())
            if mole_object.is_valid and len(mole_object.query_entry.get()) != 0:
                e = doc.createElement("NonActiveParts")
                child = doc.createElement("Query")
                child.appendChild(doc.createTextNode(query_str))
                e.appendChild(child)
                root.appendChild(e)

            e = doc.createElement("Params")
            child = doc.createElement("Cavity")

            for i in cavityParameters.keys():
                child.setAttribute(i, cavityParameters[i])

            if mole_object.ignore_het.get():
                child.setAttribute('IgnoreHETAtoms', '1')

            if mole_object.remove_hydrogens.get():
                child.setAttribute('IgnoreHydrogens', '1')

            e.appendChild(child)

            child = doc.createElement("Tunnel")

            for i in tunnelParameters.keys():
                child.setAttribute(i, tunnelParameters[i])

            child.setAttribute('WeightFunction', mole_object.weight_function.get())
            e.appendChild(child)

            root.appendChild(e)
            # endregion

            # region Export
            e = doc.createElement("Export")
            child = doc.createElement("Formats")
            child.setAttribute('PDBStructure', '0')
            child.setAttribute('PDBProfile', '1')
            child.setAttribute('CSV', '1')

            if mole_object.plugin_type == "PyMOL":
                child.setAttribute('PyMol', '1')

            else:
                child.setAttribute('Chimera', '1')

            child.setAttribute('Mesh', '0')
            child.setAttribute('ChargeSurface', '0')
            e.appendChild(child)

            child = doc.createElement("Types")
            child.setAttribute('Cavities', '0')

            if param == 'pores':
                if len(mole_object.path_starting_points_box.getvalue()) == 0 and len(mole_object.path_end_points_box.getvalue()) == 0:
                        child.setAttribute('PoresAuto', '1')

                else:
                    child.setAttribute('PoresUser', '1')

            e.appendChild(child)

            child = doc.createElement("PyMol")
            child.setAttribute('SurfaceType', 'Spheres')
            e.appendChild(child)

            #if param == 'pores':
                #e.setAttribute('Tunnels', '0')

            root.appendChild(e)


        if param == 'tunnels':
            e = doc.createElement("Origins")
            mole_object.compute_tunnels_button.config(text='Processing... Please wait.')

            if len(mole_object.start_points_box.getvalue()) == 0:
                e.setAttribute('Auto', '1')

            else:
                e.setAttribute('Auto', '0')
                child = doc.createElement("Origin")

            if len(mole_object.start_points_box.getvalue()) == 1:

                for i in mole_object.points.keys():

                    if i in mole_object.start_points_box.getvalue()[0]:
                        self.append_origin(doc, child, mole_object.points[i])
                        e.appendChild(child)

            if len(mole_object.start_points_box.getvalue()) > 1:

                for i in mole_object.points.keys():

                    for j in mole_object.start_points_box.getvalue():

                        if i in j:
                            f = doc.createElement('Origin')
                            self.append_origin(doc, f, mole_object.points[i])
                            e.appendChild(f)

                e.appendChild(child)

            root.appendChild(e)
        # endregion

        # region PathPointsBoxes
        if param == 'pores':

            if len(mole_object.path_starting_points_box.getvalue()) > 0:
                if len(mole_object.path_end_points_box.getvalue()) > 0:
                    e = doc.createElement('CustomExits')

                else:
                    mole_object.when_error('',
                                           'Only one Pore point selected. Please select at least one start pair of points'
                                           'or deselect, for automatic pores detection')
                    return

            mole_object.compute_pores_button.config(text='Processing... Please wait.')
            start_values = mole_object.path_starting_points_box.getvalue()
            stop_values = mole_object.path_end_points_box.getvalue()

            for i in start_values:  # todo here

                for j in stop_values:

                    start = mole_object.get_element(i)
                    stop = mole_object.get_element(j)

                    if start is not '42' and stop is not '42':
                        f = doc.createElement('Exit')
                        self.append_origin(doc, f, start)
                        e.appendChild(f)
                        g = doc.createElement('Exit')
                        self.append_origin(doc, g, stop)
                        e.appendChild(g)

            root.appendChild(e)

        # endregion
        with open(os.path.join(self.wd, 'params.xml'), 'w') as f:
            doc.writexml(f)

        mole_object.parent.update()
        mole_object.set_state('disabled')
        plat = str.lower(sys.platform)

        if 'win32' not in plat:

            if self.mono_test() is not 1:

                mole_object.when_error('Error',
                                       'Mono environment required for MOLE 2.0 in non-windows environment is not '
                                       'installed. Please go to www.mono-project.com and install it.')
                mole_object.rename_buttons(param)
                return

        try:
            t = threading.Thread(target=self.run_sub_process(mole_object.main_parameters[0], param))
            t.daemon = True
            t.start()
            t.join()

        except Exception as e:
            mole_object.when_error('',
                                   'An error occurred during processing. If this problem persists and you are a non-windows '
                                   'user try installing \'mono-devel\' package. If it does not help, please send the text '
                                   'below to the authors with the description \n#################\n ' + str(
                                       e))
            mole_object.rename_buttons(param)
            return

        mole_object.rename_buttons(param)
        mole_object.when_computation_done(OUTPUT, param)


    def mono_test(self):
        try:
            subprocess.Popen(['mono', '-V'])
            return 1

        except:
            return 0

    @staticmethod
    def append_origin(doc, element, point):
        for j in range(1, len(point)):
            e = doc.createElement(point[j]['Element'])
            for i in point[j].keys():
                if i != 'Element':
                    e.setAttribute(i, point[j][i])
            element.appendChild(e)
