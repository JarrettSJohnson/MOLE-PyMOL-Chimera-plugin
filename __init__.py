import MainWindow
import os
import sys


try:

    def __init__(self):
        self.menuBar.addmenuitem('Plugin', 'command',
                                 'MOLE',
                                 label='MOLE 2.5',
                                 command=lambda s=self: MainWindow.MainWindow())


    path = os.path.dirname(__file__)
    sys.path.append(path)


except (NameError,ImportError):
    pass