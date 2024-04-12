from pymatgen.ext.matproj import MPRester
from PyQt5.QtCore import QObject, pyqtSignal, QProcess
import numpy as np
from lmfit import minimize

from larch.xafs import feffrunner

structre_data_directory = "/Users/akhiltayal/Dropbox (BNL.GOV)/xfit_tools/database_directory/structure_data/"
feff_data_directory = "/Users/akhiltayal/Dropbox (BNL.GOV)/xfit_tools/database_directory/feff_data/"


API_KEY = "PvgUucO72EAOfiRI3Z2R5pr8mZAZu5ul"


class Worker_Retrive_MatProj_Data(QObject):
    finished = pyqtSignal()
    progress = pyqtSignal(int)


    def __init__(self, formula='FeO'):
        super(QObject, self).__init__()
        self.formula = formula.split(',')
        self.worker_document = {}

        self.formula_list = [f.strip() for f in self.formula]

    def run(self):
        for _formula in self.formula:

            with MPRester(API_KEY) as mpr:
                _documents = mpr.summary.search(formula=self.formula)

            for _d in _documents:
                _material_id = _d.material_id.string
                self.worker_document[_material_id] = _d
        self.finished.emit()


class Worker_Run_Feff_Calculation(QObject):
    finished = pyqtSignal()
    progress = pyqtSignal(int)

    def __init__(self, feff_path):
        super(QObject, self).__init__()
        self.feff_path = feff_path


    def run(self):
        for key in self.feff_path:
            feffrunner(folder=self.feff_path[key], feffinp='feff.inp').run()
        self.finished.emit()





class WorkerEXAFSFit(QObject):

    finished = pyqtSignal()
    progress = pyqtSignal(int)

    def __init__(self, k, data, parameters, fit_range, feff_data, kweight=2):
        super(QObject, self).__init__()
        self.k = k
        self.data = data
        self.parameters = parameters
        self.fit_range = fit_range
        self.feff_data = feff_data
        self.kweight = kweight
        self.buffer_parameters = None
        self.success = False




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


    def residual(self, params, x, data, fit_range, feff_data, k_weight=2):

        k1 = np.where(x >= fit_range[0])[0][0]
        k2 = np.where(x >= fit_range[1])[0][0]

        model = self.calculate_exafs_lmfit(x, params, feff_data, fit_range, k_weight)
        return (data[k1:k2] - model)

    def callback(self, params, iter, resid, *args, **kwargs):
        self.buffer_parameters = params

    def run(self):
        """Long-running task."""

        self.out = minimize(self.residual, params=self.parameters, method='least_squares', iter_cb=self.callback,
                        args=(self.k, self.data, self.fit_range, self.feff_data, self.kweight))

        while not self.out.success:
            print('In progress')

        self.success = True

        # self.finished.emit()
            # parametersrint()


