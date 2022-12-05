# SpinningTop
This project implements the simulation of a Lagrange spinning top (rigid, symmetric) using quaternions to compute rotations.
The spinning top is rendered in real time using OpenGL. It is possible to set various parameters at the beginning of the simulation (such as the top's initial angular velocity, mass and sizes) and to add external forces in real time. 

There are also many graphics options (to adjust brightness, display angular velocity, angular momentum and torque unitary vectors, hide or show the top's trail and all the scene components). These options are controlled through a context menu accessible with mouse right-click or through keyboard shortcuts.

The raw data produced by the simulation can be saved in .txt files and can be viewed directly in real time using the graph window, where you can select which data to be displayed and axis' range (setting it directly or zooming/traslating with mouse). All saved data are in SI units.

You should be able to build the project in Visual Studio 2019 just by downloading the repo and opening the .sln file.
Precompiled binary executable for Microsoft Windows is available in release.


To build it on Linux one of these commands should work (provided they are run in 'Spinning top' source folder):


  $ g++ main.cpp window.cpp graph.cpp `fltk-config --use-forms --use-gl  --ldstaticflags --cxxflags`
  
  or
  
  $ g++ main.cpp window.cpp graph.cpp -lGLU  -lfltk -lGL -lfltk_gl


This project relies on FLTK 1.3.6 library (www.fltk.org). FLTK' header files are included in this repository, as well as compiled x64 libraries (.lib) for Microsoft Windows.
