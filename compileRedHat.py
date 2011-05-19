#cxFreeze Compile Script
#
#
#to run type: python compile.py build
#
#

import cx_Freeze
import sys
import matplotlib

#from mvpa.suite import *

pythonLib="/usr/local/python26/lib/python2.6/site-packages/"

base = None
if sys.platform == "win32":
    base = "Win32GUI"
    includeDependencies = []
else:
    #add the dependencies bellow the exe folder by the cxFreeze
    includeDependencies = [(pythonLib+"mdp/hinet/hinet.css", "mdp/hinet/hinet.css"),
                           (pythonLib+"mdp/utils/slideshow.css", "mdp/utils/slideshow.css"),
                           ("/usr/lib64/libcblas.so","libcblas.so"),
                           ("/usr/lib64/libatlas.so","libatlas.so"),
                           ("/usr/lib64/liblapack.so","liblapack.so"),
                           ("/usr/lib64/libf77blas.so","libf77blas.so"),

#                           ("/usr/lib/libblas/libblas.so.3gf.0","libblas.so.3gf"),
#                           ("/usr/lib/lapack/liblapack.so.3gf.0","liblapack.so.3gf"),
#                           ("/usr/lib/libgfortran.so.3","libgfortran.so.3"),
                           ("/usr/lib64/libgfortran.so.1","libgfortran.so.1"),
#                           ("/etc/matplotlibrc","matplotlibrc"),
#                           ("/usr/lib/libBLT.2.4.so.8.4","libBLT.2.4.so.8.4"),
                           ("/usr/lib64/libtk8.4.so","libtk8.4.so"),
                           ("/usr/lib64/libtcl8.4.so","libtcl8.4.so"),
                           (pythonLib+"/mvpa/clfs/libsmlrc/smlrc.so","mvpa/clfs/libsmlrc/smlrc.so"),
                           #("/usr/lib64/libz.so","libz.so.1"),
                           #("/lib64/ld-2.12.1.so","ld-2.12.1.so"),
                           #("/lib64/libc.so.6","libc.so.6"),
                           #("/lib64/ld-linux-x86-64.so.2","ld-linux-x86-64.so.2"),
                           ("/lib64/libssl.so.6","libssl.so.6"),
                           ("/lib64/libcrypto.so.6","libcrypto.so.6")
                           ] #,( matplotlib.get_data_path(),"mpl-data"),]
    zipDependencies = [ (pythonLib+"/matplotlib/mpl-data/","mpl-data"),]



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
