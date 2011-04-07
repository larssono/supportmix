#cxFreeze Compile Script
#
#
#to run type: python compile.py build
#
#

import cx_Freeze
import sys
import matplotlib


base = None
if sys.platform == "win32":
    base = "Win32GUI"
    includeDependencies = []
else:
    #add the dependencies bellow the exe folder by the cxFreeze
    includeDependencies = [("/usr/share/pyshared/mdp/hinet/hinet.css", "mdp/hinet/hinet.css"),("/usr/share/pyshared/mdp/utils/slideshow.css", "mdp/utils/slideshow.css"),("/usr/lib/libblas/libblas.so.3gf.0","libblas.so.3gf"),
("/usr/lib/lapack/liblapack.so.3gf.0","liblapack.so.3gf"), ("/usr/lib/libgfortran.so.3","libgfortran.so.3"), ("/etc/matplotlibrc","matplotlibrc"),
("/usr/lib/libBLT.2.4.so.8.5","libBLT.2.4.so.8.5"), ("/usr/lib/libtk8.5.so.0","libtk8.5.so.0"), ("/usr/lib/libtcl8.5.so.0","libtcl8.5.so.0"), ("/usr/lib64/libz.so","libz.so"), ("/lib64/libc.so.6","libc.so.6"), ("/lib64/ld-linux-x86-64.so.2","ld-linux-x86-64.so.2")] #,( matplotlib.get_data_path(),"mpl-data"),]
    zipDependencies = [ ( matplotlib.get_data_path(),"mpl-data"),]



buildOptions = dict(include_files = includeDependencies,
        excludes = [""], includes = ["Tkinter","scipy","scipy.misc", "scipy.misc.common"], compressed=True, zip_includes=zipDependencies, copy_dependent_files=True)

executables = [
        cx_Freeze.Executable("SupportMix", base = base)
]

cx_Freeze.setup(
        name = "SupportMix",
        version = "0.1",
        description = "Compiled Version of Support Mix",
        executables = executables,
        options = dict(build_exe = buildOptions))
