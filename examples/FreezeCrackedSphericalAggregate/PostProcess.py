# ------------------------------------------------------------------------------
#
#  Post-processing
#
# ------------------------------------------------------------------------------

import gmsh
import os

model = gmsh.model
factory = model.geo

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)

# Let us for example include a three-dimensional scalar view:
path = os.path.dirname(os.path.abspath(__file__))
gmsh.merge('Sph2D.pos1')
gmsh.merge('Sph2D.pos2')
gmsh.merge('Sph2D.pos3')
gmsh.merge('Sph2D.pos4')
gmsh.merge('Sph2D.pos5')
gmsh.merge('Sph2D.pos6')
gmsh.merge('Sph2D.pos7')
gmsh.merge('Sph2D.pos8')
gmsh.merge('Sph2D.pos9')
gmsh.merge('Sph2D.pos10')
gmsh.merge('Sph2D.pos11')
gmsh.merge('Sph2D.pos12')
gmsh.merge('Sph2D.pos13')
gmsh.merge('Sph2D.pos14')
gmsh.merge('Sph2D.pos15')
gmsh.merge('Sph2D.pos16')
#gmsh.merge('Sph2D.pos17')
gmsh.write("iliass.pdf")
# We then set some options for the `Isosurface' plugin (which extracts an
# isosurface from a 3D scalar view), and run it:
plugin = gmsh.plugin

#plugin.setNumber("Smooth", "View", 6)
#plugin.run("Smooth")
import numpy as np


#gmsh.view.add("xxxxx")
################ MathEval (15)
plugin = gmsh.plugin
plugin.setString("MathEval", "Expression0", "v0*Cos(Atan(y/x))*Cos(Atan(y/x))+v4*Sin(Atan(y/x))*Sin(Atan(y/x))+v1*Sin(2*Atan(y/x))");
plugin.setString("MathEval", "Expression1", "(v4-v0)*Sin(Atan(y/x))*Cos(Atan(y/x))+v1*Cos(2*Atan(y/x))");
plugin.setString("MathEval", "Expression2", "v2");
plugin.setString("MathEval", "Expression3", "(v4-v0)*Sin(Atan(y/x))*Cos(Atan(y/x))+v1*Cos(2*Atan(y/x))");
plugin.setString("MathEval", "Expression4", "v0*Sin(Atan(y/x))*Sin(Atan(y/x))+v4*Cos(Atan(y/x))*Cos(Atan(y/x))-v1*Sin(2*Atan(y/x))");
plugin.setString("MathEval", "Expression5", "v5");
plugin.setString("MathEval", "Expression6", "v6");
plugin.setString("MathEval", "Expression7", "v7");
plugin.setString("MathEval", "Expression8", "v8");
plugin.setNumber("MathEval", "View", 10)
plugin.run("MathEval")
#Liquid pressure * Liquid saturation (16)
plugin = gmsh.plugin
plugin.setString("MathEval", "Expression0", "v0*w0");
plugin.setString("MathEval", "Expression1", "");
plugin.setString("MathEval", "Expression2", "");
plugin.setString("MathEval", "Expression3", "");
plugin.setString("MathEval", "Expression4", "");
plugin.setString("MathEval", "Expression5", "");
plugin.setString("MathEval", "Expression6", "");
plugin.setString("MathEval", "Expression7", "");
plugin.setString("MathEval", "Expression8", "");
plugin.setNumber("MathEval", "View", 0)
plugin.setNumber("MathEval", "OtherView", 3)
plugin.run("MathEval")
#Ice pressure * Ice saturation (17)
plugin = gmsh.plugin
plugin.setString("MathEval", "Expression0", "v0*(1-w0)");
plugin.setString("MathEval", "Expression1", "");
plugin.setString("MathEval", "Expression2", "");
plugin.setString("MathEval", "Expression3", "");
plugin.setString("MathEval", "Expression4", "");
plugin.setString("MathEval", "Expression5", "");
plugin.setString("MathEval", "Expression6", "");
plugin.setString("MathEval", "Expression7", "");
plugin.setString("MathEval", "Expression8", "");
plugin.setNumber("MathEval", "View", 6)
plugin.setNumber("MathEval", "OtherView", 3)
plugin.run("MathEval")

#Ice pressure * Ice saturation + Liquid pressure * Liquid saturation (18)
plugin = gmsh.plugin
plugin.setString("MathEval", "Expression0", "v0+w0");
plugin.setString("MathEval", "Expression1", "");
plugin.setString("MathEval", "Expression2", "");
plugin.setString("MathEval", "Expression3", "");
plugin.setString("MathEval", "Expression4", "");
plugin.setString("MathEval", "Expression5", "");
plugin.setString("MathEval", "Expression6", "");
plugin.setString("MathEval", "Expression7", "");
plugin.setString("MathEval", "Expression8", "");
plugin.setNumber("MathEval", "View", 16)
plugin.setNumber("MathEval", "OtherView", 17)
plugin.run("MathEval")

option = gmsh.option
option.setString("View[18].Name", "sl*pl+sc*pc")

# Mean stress (19)
plugin = gmsh.plugin
plugin.setString("MathEval", "Expression0", "(v0+v4+v8)/3");
plugin.setString("MathEval", "Expression1", "");
plugin.setString("MathEval", "Expression2", "");
plugin.setString("MathEval", "Expression3", "");
plugin.setString("MathEval", "Expression4", "");
plugin.setString("MathEval", "Expression5", "");
plugin.setString("MathEval", "Expression6", "");
plugin.setString("MathEval", "Expression7", "");
plugin.setString("MathEval", "Expression8", "");
plugin.setNumber("MathEval", "View", 15)
plugin.run("MathEval")

option = gmsh.option
option.setString("View[19].Name", "Mean stress")

# Effective mean stress (20)
plugin = gmsh.plugin
plugin.setString("MathEval", "Expression0", "v0+w0");
plugin.setString("MathEval", "Expression1", "");
plugin.setString("MathEval", "Expression2", "");
plugin.setString("MathEval", "Expression3", "");
plugin.setString("MathEval", "Expression4", "");
plugin.setString("MathEval", "Expression5", "");
plugin.setString("MathEval", "Expression6", "");
plugin.setString("MathEval", "Expression7", "");
plugin.setString("MathEval", "Expression8", "");
plugin.setNumber("MathEval", "View", 19)
plugin.setNumber("MathEval", "OtherView", 18)
plugin.run("MathEval")

option = gmsh.option
option.setString("View[20].Name", "Effective mean stress")

# Von-Mises stress (21)
plugin = gmsh.plugin
plugin.setString("MathEval", "Expression0", "sqrt(( (v0-v4)*(v0-v4) + (v4-v8)*(v4-v8) + (v8-v0)*(v8-v0) +6*(v1*v1+v5*v5+v6*v6)  )/2)"); #
plugin.setString("MathEval", "Expression1", "");
plugin.setString("MathEval", "Expression2", "");
plugin.setString("MathEval", "Expression3", "");
plugin.setString("MathEval", "Expression4", "");
plugin.setString("MathEval", "Expression5", "");
plugin.setString("MathEval", "Expression6", "");
plugin.setString("MathEval", "Expression7", "");
plugin.setString("MathEval", "Expression8", "");
plugin.setNumber("MathEval", "View", 10)
plugin.run("MathEval")

option = gmsh.option
option.setString("View[21].Name", "Von-Mises stress")
######################################################################
#Sigma_rr+slpl+scpc
plugin = gmsh.plugin
plugin.setString("MathEval", "Expression0", "v0+w0");
plugin.setString("MathEval", "Expression1", "");
plugin.setString("MathEval", "Expression2", "");
plugin.setString("MathEval", "Expression3", "");
plugin.setString("MathEval", "Expression4", "");
plugin.setString("MathEval", "Expression5", "");
plugin.setString("MathEval", "Expression6", "");
plugin.setString("MathEval", "Expression7", "");
plugin.setString("MathEval", "Expression8", "");
plugin.setNumber("MathEval", "View", 15)
plugin.setNumber("MathEval", "OtherView", 18)
plugin.run("MathEval")

option = gmsh.option
option.setString("View[22].Name", "Calculated cohesion")
####################################################################
#Cohesion (Sigma_rr-tan(phi)*Sigma_rt) (23)
plugin = gmsh.plugin
plugin.setString("MathEval", "Expression0", "w0-tan(35*3.14/180)*Abs(v1)");
plugin.setString("MathEval", "Expression1", "");
plugin.setString("MathEval", "Expression2", "");
plugin.setString("MathEval", "Expression3", "");
plugin.setString("MathEval", "Expression4", "");
plugin.setString("MathEval", "Expression5", "");
plugin.setString("MathEval", "Expression6", "");
plugin.setString("MathEval", "Expression7", "");
plugin.setString("MathEval", "Expression8", "");
plugin.setNumber("MathEval", "View", 15)
plugin.setNumber("MathEval", "OtherView", 18)
plugin.run("MathEval")

option = gmsh.option
option.setString("View[23].Name", "Calculated cohesion")
#####################################################################
#Cohesion (Drucker-Prager (effective mean stress vs VM stress) (24)
plugin = gmsh.plugin
plugin.setString("MathEval", "Expression0", "-0.54*v0+w0");
plugin.setString("MathEval", "Expression1", "");
plugin.setString("MathEval", "Expression2", "");
plugin.setString("MathEval", "Expression3", "");
plugin.setString("MathEval", "Expression4", "");
plugin.setString("MathEval", "Expression5", "");
plugin.setString("MathEval", "Expression6", "");
plugin.setString("MathEval", "Expression7", "");
plugin.setString("MathEval", "Expression8", "");
plugin.setNumber("MathEval", "View", 20)
plugin.setNumber("MathEval", "OtherView", 21)
plugin.run("MathEval")

option = gmsh.option
option.setString("View[24].Name", "Drucker-Prager envelope (Effective stress)")

#### cut cement paste
plugin.setString("CutParametric", "X", "0.00401*cos(u)")
plugin.setString("CutParametric", "Y", "0.00401*sin(u)")
plugin.setString("CutParametric", "Z", "0")
plugin.setNumber("CutParametric", "MinU", 1.56)
plugin.setNumber("CutParametric", "MaxU", -1.56)
plugin.setNumber("CutParametric", "NumPointsU", 40)
plugin.setNumber("CutParametric", "MinV", 0)
plugin.setNumber("CutParametric", "MaxV", 0)
plugin.setNumber("CutParametric", "NumPointsV", 1)
plugin.setNumber("CutParametric", "ConnectPoints", 0)
plugin.setNumber("CutParametric", "View", 15)
plugin.run("CutParametric")

#### cut aggregate
plugin.setString("CutParametric", "X", "0.00399*cos(u)")
plugin.setString("CutParametric", "Y", "0.00399*sin(u)")
plugin.setString("CutParametric", "Z", "0")
plugin.setNumber("CutParametric", "MinU", 1.56)
plugin.setNumber("CutParametric", "MaxU", -1.56)
plugin.setNumber("CutParametric", "NumPointsU", 40)
plugin.setNumber("CutParametric", "MinV", 0)
plugin.setNumber("CutParametric", "MaxV", 0)
plugin.setNumber("CutParametric", "NumPointsV", 1)
plugin.setNumber("CutParametric", "ConnectPoints", 0)
plugin.setNumber("CutParametric", "View", 15)
plugin.run("CutParametric")

#### cut interface
plugin.setString("CutParametric", "X", "0.004*cos(u)")
plugin.setString("CutParametric", "Y", "0.004*sin(u)")
plugin.setString("CutParametric", "Z", "0")
plugin.setNumber("CutParametric", "MinU", 1.56)
plugin.setNumber("CutParametric", "MaxU", -1.56)
plugin.setNumber("CutParametric", "NumPointsU", 40)
plugin.setNumber("CutParametric", "MinV", 0)
plugin.setNumber("CutParametric", "MaxV", 0)
plugin.setNumber("CutParametric", "NumPointsV", 1)
plugin.setNumber("CutParametric", "ConnectPoints", 0)
plugin.setNumber("CutParametric", "View", 15)
plugin.run("CutParametric")

#### cut interface
plugin.setString("CutParametric", "X", "0.004*cos(u)")
plugin.setString("CutParametric", "Y", "0.004*sin(u)")
plugin.setString("CutParametric", "Z", "0")
plugin.setNumber("CutParametric", "MinU", 1.56)
plugin.setNumber("CutParametric", "MaxU", -1.56)
plugin.setNumber("CutParametric", "NumPointsU", 400)
plugin.setNumber("CutParametric", "MinV", 0)
plugin.setNumber("CutParametric", "MaxV", 0)
plugin.setNumber("CutParametric", "NumPointsV", 1)
plugin.setNumber("CutParametric", "ConnectPoints", 0)
plugin.setNumber("CutParametric", "View", 22)
plugin.run("CutParametric")

# We finish by setting some options:

#option.setNumber("View[0].IntervalsType", 1)
#option.setNumber("View[0].NbIso", 6)
#option.setNumber("View[0].SmoothNormals", 1)
#option.setNumber("View[1].IntervalsType", 2)
#option.setNumber("View[2].IntervalsType", 2)

# show the GUI at the end
gmsh.fltk.run()

gmsh.finalize()
