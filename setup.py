import sys
import os
os.environ['TCL_LIBRARY'] = "C:\\Users\steph\\anaconda3\\tcl\\tcl8.6"
os.environ['TK_LIBRARY'] = "C:\\Users\steph\\anaconda3\\tcl\\tK8.6"
from cx_Freeze import setup, Executable
import cx_Freeze.hooks
build_exe_options = {"packages": ['pygments.lexers', 'tvtk.pyface.ui.qt4','pkg_resources._vendor',
                    'mayavi', 'traits', 'traitsui','numpy', 'matplotlib', 'math','scipy'
                    'tvtk.vtk_module','traits.api','traitsui.api','os','mayavi', 'tvtk.vtk_module',
                                  'pyface.qt','pyface.qt.QtGui','pyface.qt.QtCore','numpy'],
                     "includes":['mayavi','PyQt5'],
                    }
executables = [
    Executable('orbitalgui_June7th.py', targetName='orbitalgui_June7th.exe',base = 'Win32GUI',)
]
setup(name = "orbitalgui_June7th" ,
      version = "0.1",
      description = "",
      executables = executables,
      options={"build_exe":build_exe_options},
      )