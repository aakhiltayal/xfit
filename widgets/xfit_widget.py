# from widgets.shell_widget import Shell
from PyQt5 import uic
import os
from lmfit import Parameters
from PyQt5 import QtWidgets
import numpy as np
import pandas as pd

from PyQt5.Qt import Qt
from PyQt5.QtGui import QColor

from PyQt5.QtCore import QObject, QThread, pyqtSignal
from larch.symboltable import Group

from pyqtgraph import mkPen

from .workers import Worker_Retrive_MatProj_Data as worker_matproj
from .workers import Worker_Run_Feff_Calculation as worker_feff
from .workers import WorkerEXAFSFit as worker_fit
from .pyqtgraph_widget import create_pyqtgraph_widget
from .pyqtgraph_widget import MplCanvas

from scipy.interpolate import interp1d

from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT

from .often_used_functions import read_lineEdit_and_perform_sanity_check, get_default_feff_parameters, Shell

from pymatgen.io.feff.sets import MPXANESSet, MPEXAFSSet, FEFFDictSet

from lightshow import FEFFParameters

from xraydb import atomic_symbol

from larch.xafs import xftf

from larch.xafs.feffutils import get_feff_pathinfo

from larch.xafs.feffdat import *

ATOMIC_SYMBOL_DICT = {'element': [atomic_symbol(i) for i in range(20, 93)]}
EDGES = {'edges': ['K', 'L3', 'L2', 'L1', 'M5']}


STRUCTURE_FOLDER = "/Users/akhiltayal/Dropbox (BNL.GOV)/xfit_tools/database_directory/structure_data/"


PATH = '/Users/akhiltayal/Dropbox (BNL.GOV)/xfit_tools/ui/'

class plot_app(*uic.loadUiType(PATH + 'plot_app.ui')):
    def __init__(self):
        super().__init__()
        self.setupUi(self)
        self.pushButton_search.clicked.connect(self.search_the_structure_from_materials_project)
        self.material_project_data = {}
        self.documents = None
        self.create_plot_item_for_simulated_chi()
        self.create_plot_item_for_simulated_ft()
        self.create_plot_item_for_fit_chi()
        self.create_plot_item_for_fit_ft()
        self.pushButton_run_feff.clicked.connect(self.perform_feff_calculation)
        self.pushButton_list_path.clicked.connect(self.populate_scattering_paths)
        self.pushButton_add_scattering_path.clicked.connect(self.plot_scattering_path)
        self.pushButton_add_shells_to_fit.clicked.connect(self.add_selected_paths_tofit)
        self.pushButton_load_chi.clicked.connect(self.load_chi_data)
        self.pushButton_make_ft.clicked.connect(self.plot_raw_chi_ft)
        self.pushButton_fit.clicked.connect(self.perform_exafs_fit)
        self.pushButton_save_current_session.clicked.connect(self.create_hdf5_file_current_session)
        self.structure_found = False
        self.path_added = False
        self.shells = {}
        self.parameter = Parameters()
        self.treeWidget_paths_dict = {}

        windows = ['sine', 'hanning', 'parzen', 'welch', 'gaussian', 'kaiser']

        for key in ["", '_sim']:
            getattr(self, 'comboBox_window' + key).addItems(windows)

    def search_the_structure_from_materials_project(self):
        _formula = self.lineEdit_formula.text()
        self.search_materials_structure(_formula)


    def search_materials_structure(self, formula='FeO'):
        print("Searching Material...\n")
        self.thread_worker_matproj = QThread()
        self.worker_matproj = worker_matproj(formula=formula)
        self.worker_matproj.moveToThread(self.thread_worker_matproj)
        self.thread_worker_matproj.started.connect(self.worker_matproj.run)
        self.worker_matproj.finished.connect(self.thread_worker_matproj.quit)
        self.worker_matproj.finished.connect(self.populate_materials_structure)
        self.thread_worker_matproj.finished.connect(self.get_finished_status)
        self.worker_matproj.finished.connect(self.worker_matproj.deleteLater)
        self.thread_worker_matproj.start()

    def populate_materials_structure(self):
        if self.worker_matproj.worker_document is not None:
            self.documents = self.worker_matproj.worker_document

        if len(self.documents) > 0:

            _labels = ['mp-ID', 'Formula', 'structure', 'E full(eV)', "Experimentally Observed"]
            self.clear_treeWidget(tree_widget=self.treeWidget_structure, labels=_labels)
            _parent = self.treeWidget_structure
            self._treeWidget = {}
            print("Loading Structure...\n")
            for key, doc in self.documents.items():
                _name_list = [doc.material_id.string,
                              doc.formula_pretty,
                              f"{doc.symmetry.crystal_system}",
                              f"{doc.energy_above_hull:2.3f}",
                              f"{not doc.theoretical}"]
                _status, _status_experimentally_observed = self.check_if_feff_files_will_be_good(doc) #Must be removed when pymatgen update the code
                self._treeWidget[doc.material_id.string] = self._make_item(parent=_parent, item_list=_name_list)
                if _status:
                    self._treeWidget[doc.material_id.string].setBackground(0, QColor('yellow'))
                if _status_experimentally_observed:
                    self._treeWidget[doc.material_id.string].setBackground(4, QColor('lime'))
            self.structure_found = True
        else:
            print('No structure found\n')

    def check_if_feff_files_will_be_good(self, document):
        for el, amt in document.composition.items():
            if amt == 1:
                status = False
            else:
                status = True
        return status, not document.theoretical

    def clear_treeWidget(self, tree_widget=None, labels=None):
        tree_widget.clear()
        if labels is not None:
            tree_widget.setHeaderLabels(labels)
            for i in range(len(labels)):
                tree_widget.setColumnWidth(i, 150)
            tree_widget.setSortingEnabled(True)

    def _make_item(self, parent=None, item_list=None):
        _item = QtWidgets.QTreeWidgetItem(parent, item_list)
        _item.setFlags(_item.flags() | Qt.ItemIsTristate | Qt.ItemIsUserCheckable)
        _item.setCheckState(0, Qt.Unchecked)
        return _item

    def get_finished_status(self):
        print('Search complete\n')

    def plot_scattering_path(self):

        self.canvas_sim_chi.axes.cla()
        self.canvas_sim_ft.axes.cla()

        self.treeWidget_paths_dict['selected_feff_data'] = {}
        __kweight = self.spinBox_kweight_sim.value()
        __window = self.comboBox_window_sim.currentText()
        __kmin = self.doubleSpinBox_kmin_sim.value()
        __kmax = self.doubleSpinBox_kmax_sim.value()

        __sum_chi = 0
        __sum_k = 0
        count = 0
        for key in self.treeWidget_paths_dict['widgets'].keys():
            self.treeWidget_paths_dict['selected_feff_data'][key] = {}
            for k, value in self.treeWidget_paths_dict['widgets'][key].items():
                if value.checkState(0) == 2:
                    __path = os.path.join(self.treeWidget_paths_dict['paths'][key], k)
                    _feffpath = feffpath(filename=__path)
                    _group = ff2chi([_feffpath])
                    __sum_k += _group.k
                    __sum_chi += _group.chi
                    count += 1
                    xftf(_group.k, _group.chi, kweight=__kweight, window=__window, kmin=__kmin, kmax=__kmax, group=_group)
                    self.treeWidget_paths_dict['selected_feff_data'][key][k] = [_feffpath, _group]
                    self.canvas_sim_chi.axes.plot(_group.k, _group.chi*_group.k**__kweight, label=f"{key} {k}")
                    self.canvas_sim_ft.axes.plot(_group.r, _group.chir_mag, label=f"{key} {k}")


        __sum = Group()
        __sum.k = __sum_k/count
        __sum.chi = __sum_chi
        xftf(__sum.k, __sum.chi, kweight=__kweight, window=__window, kmin=__kmin, kmax=__kmax, group=__sum)
        self.canvas_sim_chi.axes.plot(__sum.k, __sum.chi * __sum.k ** __kweight, label=f"sum", linewidth=2)
        self.canvas_sim_ft.axes.plot(__sum.r, __sum.chir_mag, label=f"sum", linewidth=2)

        self.canvas_sim_chi.axes.legend()
        self.canvas_sim_chi.draw()
        self.canvas_sim_ft.axes.legend()
        self.canvas_sim_ft.draw()

    def create_fourier_trasform(self, ):
        pass








    def create_chi_and_ft(self):
        self.sim_chi_ft = {}
        for dic in self.feff_paths.keys():
            for key in self.feff_paths[dic].keys():
                self.sim_chi_ft[dic][key] = ff2chi(self.feff_paths[dic][key])

            self.sim_chi_ft[dic]['summation'] = ff2chi(self.feff_paths[dic])

        _feff = []
        for key in self.feff_paths.keys():
            _feff.append(self.feff_paths[key].values())
        self.sim_chi_ft['summation'] = ff2chi(_feff)








    def make_feff_folder(self, absorbing_atom_valid_input, edge_valid_input, absorbing_atom, edge):
        task_finished = False
        if absorbing_atom_valid_input and edge_valid_input:
            self.treeWidget_paths_dict['paths'] = {}
            for key in self._treeWidget.keys():
                if self._treeWidget[key].checkState(0) == 2:
                    self.treeWidget_paths_dict['paths'][key] = os.path.join(STRUCTURE_FOLDER, self.documents[key].formula_pretty, key)
                    default_feff_parameters = get_default_feff_parameters(title=self.documents[key].material_id.string)
                    feff_dict = FEFFDictSet(absorbing_atom=absorbing_atom,
                                            structure=self.documents[key].structure,
                                            spectrum="EXAFS",
                                            config_dict=default_feff_parameters,
                                            edge=edge,
                                            radius=7.5)
                    feff_dict.write_input(self.treeWidget_paths_dict['paths'][key])
                    task_finished = True
        if not task_finished:
            print("Select at least one Structure!")


    def run_feff(self, path_dict):
        print("Searching Material...\n")

        self.thread_feff_runner = QThread()
        self.worker_feff = worker_feff(feff_path=path_dict)
        self.worker_feff.moveToThread(self.thread_feff_runner)
        self.thread_feff_runner.started.connect(self.worker_feff.run)
        self.worker_feff.finished.connect(self.thread_feff_runner.quit)
        # self.worker_feff.finished.connect(self.populate_materials_structure)
        # self.thread_feff_runner.finished.connect(self.get_finished_status)
        self.worker_feff.finished.connect(self.worker_feff.deleteLater)

        self.thread_feff_runner.start()


    def populate_scattering_paths(self):

        _labels = ['ID', 'Feff File', 'R', 'Coordination Number', 'Weight', 'Legs', 'Geometry']
        self.clear_treeWidget(tree_widget=self.treeWidget_scattering_path, labels=_labels)
        _parent = self.treeWidget_scattering_path

        self.treeWidget_paths_dict['widgets'] = {}
        self.treeWidget_paths_dict['feff_info'] = {}

        for key in self.treeWidget_paths_dict['paths'].keys():
            self.treeWidget_paths_dict['feff_info'][key] = get_feff_pathinfo(self.treeWidget_paths_dict['paths'][key])

            self.treeWidget_paths_dict['widgets'][key] = {}
            for path in self.treeWidget_paths_dict['feff_info'][key].paths:

                _name_list = [f"{key}",
                              f"{path.filename}",
                              f"{path.reff:1.3f}",
                              f"{path.degen:.0f}",
                              f"{path.cwratio:3.2f}",
                              f"{path.nleg}",
                              f"{path.geom}"]

                self.treeWidget_paths_dict['widgets'][key][path.filename] = self._make_item(parent=_parent, item_list=_name_list)
        self.path_added = True



    def perform_feff_calculation(self):
        absorbing_atom, absorbing_atom_valid_input = read_lineEdit_and_perform_sanity_check(lineEdit=self.lineEdit_absorbing_atom,
                                                                                            reference_dictionary=ATOMIC_SYMBOL_DICT,
                                                                                            message="Provide Valid Absorbing atom!")
        edge, edge_valid_input = read_lineEdit_and_perform_sanity_check(lineEdit=self.lineEdit_edge,
                                                                        reference_dictionary=EDGES,
                                                                        message="Provide Valid Edge!")

        self.make_feff_folder(absorbing_atom_valid_input=absorbing_atom_valid_input,
                              edge_valid_input=edge_valid_input,
                              absorbing_atom=absorbing_atom,
                              edge=edge)

        self.run_feff(self.treeWidget_paths_dict['paths'])





    def create_plot_item_for_simulated_chi(self, number_of_references=1):

        self.canvas_sim_chi = MplCanvas(parent=self)
        _navi = NavigationToolbar2QT(self.canvas_sim_chi, parent=self)
        self.verticalLayout_sim_chi.addWidget(_navi)
        self.verticalLayout_sim_chi.addWidget(self.canvas_sim_chi)

    def create_plot_item_for_simulated_ft(self):

        self.canvas_sim_ft = MplCanvas(parent=self)
        _navi = NavigationToolbar2QT(self.canvas_sim_ft, parent=self)
        self.verticalLayout_sim_ft.addWidget(_navi)
        self.verticalLayout_sim_ft.addWidget(self.canvas_sim_ft)


    def create_plot_item_for_fit_chi(self):
        plt_item, _ref = create_pyqtgraph_widget(layout=self.verticalLayout_chi,
                                                 title="EXAFS Chi",
                                                 number_of_references=2)

        self.raw_chi_ref = _ref[1]
        self.fit_chi_ref = _ref[2]

    def create_plot_item_for_fit_ft(self):

        plt_item, _ref = create_pyqtgraph_widget(layout=self.verticalLayout_ft,
                                                 title="Fourier Transform",
                                                 number_of_references=4)

        self.raw_ft_mag_ref = _ref[1]
        self.raw_ft_img_ref = _ref[2]
        self.fit_ft_mag_ref = _ref[3]
        self.fit_ft_img_ref = _ref[4]






    def load_chi_data(self):
        _filespath = QtWidgets.QFileDialog.getOpenFileNames(parent=self,
                                                           directory="/Users/akhiltayal/Library/CloudStorage/GoogleDrive-akhil.tayal.bnl@gmail.com/My Drive/Demo/",
                                                           options=QtWidgets.QFileDialog.DontUseNativeDialog)[0][0]

        print(_filespath)

        col_names = ['k', 'chi', 'chik', 'chik2', 'chik3', 'win', 'energy']
        self.chi_data = pd.read_csv(_filespath, comment='#', sep='\s+', names=col_names)
        self.plot_raw_chi_ft()




    def default_feff_parameters_for_exafs(self):
        self.feff_params = FEFFParameters(
            cards={
                "S02": "1",
                "EDGE": "K",
                "CONTROL": "1 1 1 1 1 1",
                "PRINT": "1 0 0 0 0 3",
                "RPATH": "6.0",
                "EXAFS": "20"
            },
            edge="K",
            radius=7.5)
        return self.feff_params

    def add_selected_paths_tofit(self):

        self.clear_shell_widgets()
        self.shells = {}
        count = 0
        for key in self.treeWidget_paths_dict['selected_feff_data'].keys():
            for i, k in enumerate(self.treeWidget_paths_dict['selected_feff_data'][key].keys()):
                self.shells['Shell_' + str(count+1)] = {}
                self.shells['Shell_' + str(count+1)]['widget'] = Shell()
                self.horizontalLayout_param.insertWidget(i, self.shells['Shell_' + str(count+1)]['widget'])
                self.shells['Shell_' + str(count+1)]['widget'].groupBox.setTitle(f'Shell {count+1} {key} {k} ')
                self.shells['Shell_' + str(count+1)]['parameter'] = self.treeWidget_paths_dict['selected_feff_data'][key][k][0]
                count += 1

        self.spinBox_shells.setValue(count)
        self.populate_shells_with_default_params(self.shells.keys())
        self.read_feff_and_interpolate(self.shells.keys())

    def populate_shells_with_default_params(self, keys):

        for i, key in enumerate(keys):
            if i == 0:
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + 'r_value').setValue(self.shells[key]['parameter'].reff)
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + 'r_max').setValue(self.shells[key]['parameter'].reff*1.05)
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + 'r_min').setValue(self.shells[key]['parameter'].reff*0.95)
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + 'n_value').setValue(self.shells[key]['parameter'].degen)
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + 'n_max').setValue(self.shells[key]['parameter'].degen*1.05)
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + 'n_min').setValue(self.shells[key]['parameter'].degen*0.95)
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + 'ss_value').setValue(0.05)
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + 'ss_max').setValue(0.2)
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + 'ss_min').setValue(0.01)
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + 'c3_value').setValue(0)
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + 'c3_max').setValue(0.02)
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + 'c3_min').setValue(0)
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + 'c4_value').setValue(0)
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + 'c4_max').setValue(0.02)
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + 'c4_min').setValue(0)
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + 'e_value').setValue(0)
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + 'e_max').setValue(20)
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + 'e_min').setValue(-20)
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + 's02_value').setValue(0.8)
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + 's02_max').setValue(1.0)
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + 's02_min').setValue(0.7)
            else:
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + 'r_value').setValue(self.shells[key]['parameter'].reff)
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + 'r_max').setValue(self.shells[key]['parameter'].reff * 1.05)
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + 'r_min').setValue(self.shells[key]['parameter'].reff * 0.95)
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + 'n_value').setValue(self.shells[key]['parameter'].degen)
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + 'n_max').setValue(self.shells[key]['parameter'].degen * 1.05)
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + 'n_min').setValue(self.shells[key]['parameter'].degen * 0.95)
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + 'ss_value').setValue(0.05)
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + 'ss_max').setValue(0.2)
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + 'ss_min').setValue(0.01)
                getattr(self.shells[key]['widget'], 'lineEdit_' + 'ss').setText(f'ss_{i-1:d}')
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + 'c3_value').setValue(0)
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + 'c3_max').setValue(0.02)
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + 'c3_min').setValue(0)
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + 'c4_value').setValue(0)
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + 'c4_max').setValue(0.02)
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + 'c4_min').setValue(0)
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + 'e_value').setValue(0)
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + 'e_max').setValue(20)
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + 'e_min').setValue(-20)
                getattr(self.shells[key]['widget'], 'lineEdit_' + 'e').setText(f'e_0')
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + 's02_value').setValue(0.8)
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + 's02_max').setValue(1.0)
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + 's02_min').setValue(0.7)
                getattr(self.shells[key]['widget'], 'lineEdit_' + 's02').setText(f's02_0')

    def clear_shell_widgets(self):
        if self.horizontalLayout_param.count() > 1:
            index = self.horizontalLayout_param.count() - 2
            while(index >= 0):
                widget = self.horizontalLayout_param.itemAt(index).widget()
                widget.setParent(None)
                index -= 1



    def perform_exafs_fit(self):
        self.pushButton_fit.setEnabled(False)
        self.params = Parameters()
        self.params.clear()

        # self.shells_feff_data = {}


        for i, key in enumerate(self.shells.keys()):

            for par in ['r', 'n', 'ss', 'c3', 'c4', 'e', 's02']:
                _value = getattr(self.shells[key]['widget'], 'doubleSpinBox_' + par + '_value').value()
                _min = getattr(self.shells[key]['widget'], 'doubleSpinBox_' + par + '_min').value()
                _max = getattr(self.shells[key]['widget'], 'doubleSpinBox_' + par + '_max').value()
                _vary = getattr(self.shells[key]['widget'], 'checkBox_' + par).isChecked()
                _expr = getattr(self.shells[key]['widget'], 'lineEdit_' + par).text()

                if i == 0 or _expr == "":
                    _expr = None

                self.params.add(f'{par}_{i:d}', value=_value, min=_min, max=_max, vary=_vary, expr=_expr)

        kmin = self.doubleSpinBox_kmin.value()
        kmax = self.doubleSpinBox_kmax.value()

        self.kweight = self.spinBox_kweight.value()

        # for key in self.shells.keys():
        #     self.shells_feff_data[key] = self.feff_data[key]



        self.fit_range = [kmin, kmax]

        self.fitting_in_progress = True



        self.shells_feff_data = [self.shells[key]['feff_data'] for key in self.shells.keys()]

        self.thread_fit = QThread()
        self.worker_fit = worker_fit(k=self.chi_data['k'],
                                     data=self.chi_data['for_fit'],
                                     parameters=self.params,
                                     fit_range=self.fit_range,
                                     feff_data=self.shells_feff_data,
                                     kweight=self.kweight)

        self.worker_fit.run()

        while not self.worker_fit.success:
            print('Fitting in progress')

        if self.worker_fit.success:
            print('Fitting Done')
            self.plot_fit_chi_ft(parameters=self.worker_fit.out.params)
            self.populate_shells_with_fit_params(parameters=self.worker_fit.out.params)


        # self.worker_fit.moveToThread(self.thread_fit)
        # self.thread_fit.started.connect(self.worker_fit.run)
        # self.worker_fit.finished.connect(self.thread_fit.quit)
        # self.worker_fit.finished.connect(self.worker_fit.deleteLater)
        # self.thread_fit.finished.connect(self.get_finished_status)
        # self.thread_fit.finished.connect(self.thread_fit.deleteLater)
        # # self.worker.progress.connect(self.reportProgress)
        # self.thread_fit.start()

        print(self.fitting_in_progress)

    def calculate_k(self, k, delE0):
        k_new = np.sqrt(k ** 2 - 0.26246592 * delE0)
        return k_new

    def calculate_exafs_lmfit(self, k, parameters, feff_data, fit_range, k_weight=2):

        parameter_dict = parameters.valuesdict()

        if fit_range is None:
            k1 = 0
            k2 = len(k)
        else:
            k1 = np.where(k >= fit_range[0])[0][0]
            k2 = np.where(k >= fit_range[1])[0][0]

        exafs = 0
        for i, data in enumerate(feff_data):

            feff = data
            R = parameter_dict['r_' + str(i)]
            N = parameter_dict['n_' + str(i)]
            sigma = parameter_dict['ss_' + str(i)]
            delE0 = parameter_dict['e_' + str(i)]
            C3 = parameter_dict['c3_' + str(i)]
            C4 = parameter_dict['c4_' + str(i)]
            S02 = parameter_dict['s02_' + str(i)]

            k_new = self.calculate_k(k, delE0)
            exafs = exafs + (S02 * N *
                             feff.mag_feff.iloc[k1:k2] / (k_new[k1:k2] * R ** 2) *
                             np.sin(2 * k_new[k1:k2] * R - ((4 / 3) * C3 * k_new[k1:k2] ** 3) + feff.pha_feff.iloc[
                                                                                                k1:k2] + feff.real_phc.iloc[
                                                                                                         k1:k2]) *
                             np.exp(-2 * k[k1:k2] ** 2 * sigma ** 2) *
                             np.exp((2 / 3) * k[k1:k2] ** 4 * C4) *
                             np.exp(-2 * R / feff.lam.iloc[k1:k2])) * (k_new[k1:k2] ** k_weight)
        return exafs


    def get_finished_status(self):
        print('Fitting finished')
        self.fitting_in_progress = False
        self.pushButton_fit.setEnabled(True)



    def plot_fit_chi_ft(self, parameters=None):

        _parameters = parameters

        exafs_fit = self.calculate_exafs_lmfit(self.chi_data['k'],
                                               _parameters,
                                               self.shells_feff_data,
                                               fit_range=None,
                                               k_weight=self.kweight)
        _buffer = Group()

        xftf(self.chi_data['k'], np.nan_to_num(exafs_fit), group=_buffer, kweight=0, window='sine', kmin=3, kmax=14)

        pen = mkPen(color='r', width=2, style=Qt.PenStyle.SolidLine)


        self.fit_chi_ref.setData(self.chi_data['k'], np.nan_to_num(exafs_fit), pen=pen)

        self.fit_ft_mag_ref.setData(_buffer.r, _buffer.chir_mag, pen=pen)
        self.fit_ft_img_ref.setData(_buffer.r, _buffer.chir_im, pen=pen)


    #
    #
    # def update_plot(self):
    #     if self.fitting_in_progress and (self.worker.buffer_parameters is not None):
    #         self.plot_fit_chi_ft(parameters=self.worker.buffer_parameters)
    #         self.populate_shells_with_fit_params(parameters=self.worker.buffer_parameters)
    #
    #
    #
    def populate_shells_with_fit_params(self, parameters=None):

        _parameters = parameters

        _variables = ['r', 'n', 'ss', 'c3', 'c4', 'e', 's02']

        for i, key in enumerate(self.shells.keys()):
            for _var in _variables:
                getattr(self.shells[key]['widget'], 'doubleSpinBox_' + _var + '_value').setValue(_parameters.valuesdict()[f'{_var}_{i:d}'])

        self.pushButton_fit.setEnabled(True)


    def add_and_remove_shells(self):
        pass


    def read_feff_and_interpolate(self, shell_keys):
        new_k = np.arange(0, 20.01, 0.05)
        keys_data = ['mag_feff', 'pha_feff', 'real_phc', 'lam']

        _buffer = {}
        _buffer['k'] = new_k
        for shell_key in shell_keys:
            for key_data in keys_data:
                x = getattr(self.shells[shell_key]['parameter']._feffdat, 'k')
                y = getattr(self.shells[shell_key]['parameter']._feffdat, key_data)
                function = interp1d(x, y, kind='cubic')
                _buffer[key_data] = function(new_k)
            self.shells[shell_key]['feff_data'] = pd.DataFrame(_buffer)



    def plot_raw_chi_ft(self):

        __kweight = self.spinBox_kweight.value()
        __window = self.comboBox_window.currentText()
        __kmin = self.doubleSpinBox_kmin.value()
        __kmax = self.doubleSpinBox_kmax.value()

        __chi = self.chi_data['chi'] *self.chi_data['k']**__kweight
        self.chi_data['for_fit'] = __chi
        __group = Group()
        xftf(self.chi_data['k'], self.chi_data['chi'], kweight=__kweight, kmin=__kmin, kmax=__kmax, window=__window, group=__group)

        pen = mkPen(color='b', width=2, style=Qt.PenStyle.SolidLine)

        self.raw_chi_ref.setData(self.chi_data['k'], __chi, pen=pen, clear=True)
        self.raw_ft_mag_ref.setData(__group.r, __group.chir_mag, pen=pen, clear=True)
        self.raw_ft_img_ref.setData(__group.r, __group.chir_im, pen=pen, clear=True)

        self.fit_chi_ref.clear()
        self.fit_ft_mag_ref.clear()
        self.fit_ft_img_ref.clear()



    def create_hdf5_file_current_session(self):
        pass









