import sys
from PyQt5 import QtWidgets, QtCore

sys.path.append('/Users/akhiltayal/Dropbox (BNL.GOV)/xfit_tools/')

from widgets.xfit_widget import plot_app


app = QtWidgets.QApplication([sys.argv])

window = plot_app()

def xfit():
    window.show()


xfit()



# app = QtWidgets.QApplication([])
# app.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)

# xfit2_gui = plot_app()

# if __name__ == '__main__':
#     app = QtWidgets.QApplication([])
#     window = plot_app()
#     window.show()
#     sys.exit(app.exec_())


# def xfit():
#     print('test')
#     xfit2_gui.show()
#     app.exec_()
#
# xfit()


