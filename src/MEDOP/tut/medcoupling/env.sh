# This configuration comes in extension with the standard configuration
# defined by the file envlight.sh for image processing and math
# processing.

#------ Python Imaging Library ------
PILDIR=/opt/programs/salome/workspace/prerequisites/install/Pil-117-py266-tcl859-tk859
export PATH=${PILDIR}/bin:$PATH
export PYTHONPATH=${PILDIR}/lib/python2.6/site-packages:${PYTHONPATH}

#------ Numpy ------
NUMPY_DIR=/opt/programs/salome/workspace/prerequisites/install/Numpy-151-py266-la331
export PATH=${NUMPY_DIR}/bin:${PATH}
export PYTHONPATH=${NUMPY_DIR}/lib/python2.6/site-packages:${PYTHONPATH}

# ------ Scipy ------
SCIPY_DIR=/opt/programs/salome/workspace/prerequisites/install/Scipy-090-py266-la331-sw204-nu151
export PYTHONPATH=${SCIPY_DIR}/lib/python2.6/site-packages:${PYTHONPATH}

# WARN: Matplot, sip and pyqt are requires for the plotter used in
# lagrange.py (could be optional)

# ------ Matplot ----
MATPLOT_DIR=/opt/programs/salome/workspace/prerequisites/install/Matplotlib-110-py266-set06c11-num151
export PYTHONPATH=${MATPLOT_DIR}/lib/python2.6/site-packages:${PYTHONPATH}

#------ sip ------
SIPDIR=/opt/programs/salome/workspace/prerequisites/install/Sip-4132-py266
export PATH=${SIPDIR}/bin:${PATH}
export PYTHONPATH=${SIPDIR}/lib/python2.6/site-packages:${PYTHONPATH}
export LD_LIBRARY_PATH=${SIPDIR}/lib/python2.6/site-packages:${LD_LIBRARY_PATH}

PYQTDIR=/opt/programs/salome/workspace/prerequisites/install/Pyqt-491p1-py266-qt463p2-sip4132
#export PYQT_SIPS=${PYQTDIR}/share/sip
#export PYUIC=${PYQTDIR}/bin/pyuic4
export PYTHONPATH=${PYQTDIR}/lib/python2.6/site-packages:${PYTHONPATH}
export PATH=${PYQTDIR}/bin:${PATH}
