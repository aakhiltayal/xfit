import pyqtgraph as pg
import numpy as np

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure


pg.setConfigOption('leftButtonPan', False)
pg.setConfigOption('background', 'white')
pg.setConfigOption('foreground', 'k')


class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None):
        super(FigureCanvasQTAgg, self).__init__()
        fig = Figure()
        self.axes = fig.add_subplot(111)
        super().__init__(fig)



def create_pyqtgraph_widget(layout=None,
                            title="pyqtgraph widget",
                            number_of_references=1):
    window = pg.GraphicsLayoutWidget()
    layout.addWidget(window)
    plot_item = window.addPlot()
    plot_item.setTitle(title, size='18', color='k')
    plot_item.addLegend()

    _ref = np.arange(1, number_of_references+1)

    references = {}
    for i in _ref:
        references[i] = plot_item.plot()

    return plot_item, references


# def create_chi_plot(self):
#
#     self.win_chi = pg.GraphicsLayoutWidget()
#     self.verticalLayout_chi.addWidget(self.win_chi)
#     self.plot_chi = self.win_chi.addPlot()
#     self.plot_chi.setTitle('Chi plot')
#
#     self.ref_raw_chi = self.plot_chi.plot()
#     self.ref_fit_chi = self.plot_chi.plot()
#
#
# def create_ft_plot(self):
#     self.win_ft = pg.GraphicsLayoutWidget()
#     self.verticalLayout_ft.addWidget(self.win_ft)
#     self.plot_ft = self.win_ft.addPlot()
#     self.plot_ft.setTitle('FT plot')
#
#     self.ref_raw_ft_mag = self.plot_ft.plot()
#     self.ref_raw_ft_img = self.plot_ft.plot()
#
#     self.ref_fit_ft_mag = self.plot_ft.plot()
#     self.ref_fit_ft_img = self.plot_ft.plot()