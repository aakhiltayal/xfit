from PyQt5 import uic

path = '/Users/akhiltayal/Dropbox (BNL.GOV)/xfit_tools/ui/'

class Shell(*uic.loadUiType(path + 'shell2.ui')):
    def __init__(self):
        super().__init__()
        self.setupUi(self)


# if __name__ == '__main__':
#     app = QtWidgets.QApplication([])
#     window = Shell()
#     window.show()
#     exit(app.exec_())