# -*- coding: mbcs -*-
# Do not delete the following import lines
from abaqus import *
from abaqusConstants import *
import __main__

def Macro_postprocess():
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    odb = session.odbs['C:/Users/aero_waas_2018_1/Documents/GitHub/composite-multiscale/micro_ruc/prepare_RUC_model/postprocessing/RUC_multi_strain_original.odb']
    session.XYDataFromHistory(name='RF1 PI: rootAssembly N: 1 NSET SETRPX-1', 
        odb=odb, 
        outputVariableName='Reaction force: RF1 PI: rootAssembly Node 1 in NSET SETRPX', 
        steps=('Step-1', ), __linkedVpName__='Viewport: 1')
    session.XYDataFromHistory(name='RF1 PI: rootAssembly N: 2 NSET SETRPY-1', 
        odb=odb, 
        outputVariableName='Reaction force: RF1 PI: rootAssembly Node 2 in NSET SETRPY', 
        steps=('Step-1', ), __linkedVpName__='Viewport: 1')
    session.XYDataFromHistory(name='RF1 PI: rootAssembly N: 3 NSET SETRPZ-1', 
        odb=odb, 
        outputVariableName='Reaction force: RF1 PI: rootAssembly Node 3 in NSET SETRPZ', 
        steps=('Step-1', ), __linkedVpName__='Viewport: 1')
    session.XYDataFromHistory(name='RF2 PI: rootAssembly N: 1 NSET SETRPX-1', 
        odb=odb, 
        outputVariableName='Reaction force: RF2 PI: rootAssembly Node 1 in NSET SETRPX', 
        steps=('Step-1', ), __linkedVpName__='Viewport: 1')
    session.XYDataFromHistory(name='RF2 PI: rootAssembly N: 2 NSET SETRPY-1', 
        odb=odb, 
        outputVariableName='Reaction force: RF2 PI: rootAssembly Node 2 in NSET SETRPY', 
        steps=('Step-1', ), __linkedVpName__='Viewport: 1')
    session.XYDataFromHistory(name='RF2 PI: rootAssembly N: 3 NSET SETRPZ-1', 
        odb=odb, 
        outputVariableName='Reaction force: RF2 PI: rootAssembly Node 3 in NSET SETRPZ', 
        steps=('Step-1', ), __linkedVpName__='Viewport: 1')
    session.XYDataFromHistory(name='RF3 PI: rootAssembly N: 1 NSET SETRPX-1', 
        odb=odb, 
        outputVariableName='Reaction force: RF3 PI: rootAssembly Node 1 in NSET SETRPX', 
        steps=('Step-1', ), __linkedVpName__='Viewport: 1')
    session.XYDataFromHistory(name='RF3 PI: rootAssembly N: 2 NSET SETRPY-1', 
        odb=odb, 
        outputVariableName='Reaction force: RF3 PI: rootAssembly Node 2 in NSET SETRPY', 
        steps=('Step-1', ), __linkedVpName__='Viewport: 1')
    session.XYDataFromHistory(name='RF3 PI: rootAssembly N: 3 NSET SETRPZ-1', 
        odb=odb, 
        outputVariableName='Reaction force: RF3 PI: rootAssembly Node 3 in NSET SETRPZ', 
        steps=('Step-1', ), __linkedVpName__='Viewport: 1')
    session.XYDataFromHistory(name='U1 PI: rootAssembly N: 1 NSET SETRPX-1', 
        odb=odb, 
        outputVariableName='Spatial displacement: U1 PI: rootAssembly Node 1 in NSET SETRPX', 
        steps=('Step-1', ), __linkedVpName__='Viewport: 1')
    session.XYDataFromHistory(name='U1 PI: rootAssembly N: 2 NSET SETRPY-1', 
        odb=odb, 
        outputVariableName='Spatial displacement: U1 PI: rootAssembly Node 2 in NSET SETRPY', 
        steps=('Step-1', ), __linkedVpName__='Viewport: 1')
    session.XYDataFromHistory(name='U1 PI: rootAssembly N: 3 NSET SETRPZ-1', 
        odb=odb, 
        outputVariableName='Spatial displacement: U1 PI: rootAssembly Node 3 in NSET SETRPZ', 
        steps=('Step-1', ), __linkedVpName__='Viewport: 1')
    session.XYDataFromHistory(name='U2 PI: rootAssembly N: 1 NSET SETRPX-1', 
        odb=odb, 
        outputVariableName='Spatial displacement: U2 PI: rootAssembly Node 1 in NSET SETRPX', 
        steps=('Step-1', ), __linkedVpName__='Viewport: 1')
    session.XYDataFromHistory(name='U2 PI: rootAssembly N: 2 NSET SETRPY-1', 
        odb=odb, 
        outputVariableName='Spatial displacement: U2 PI: rootAssembly Node 2 in NSET SETRPY', 
        steps=('Step-1', ), __linkedVpName__='Viewport: 1')
    session.XYDataFromHistory(name='U2 PI: rootAssembly N: 3 NSET SETRPZ-1', 
        odb=odb, 
        outputVariableName='Spatial displacement: U2 PI: rootAssembly Node 3 in NSET SETRPZ', 
        steps=('Step-1', ), __linkedVpName__='Viewport: 1')
    session.XYDataFromHistory(name='U3 PI: rootAssembly N: 1 NSET SETRPX-1', 
        odb=odb, 
        outputVariableName='Spatial displacement: U3 PI: rootAssembly Node 1 in NSET SETRPX', 
        steps=('Step-1', ), __linkedVpName__='Viewport: 1')
    session.XYDataFromHistory(name='U3 PI: rootAssembly N: 2 NSET SETRPY-1', 
        odb=odb, 
        outputVariableName='Spatial displacement: U3 PI: rootAssembly Node 2 in NSET SETRPY', 
        steps=('Step-1', ), __linkedVpName__='Viewport: 1')
    session.XYDataFromHistory(name='U3 PI: rootAssembly N: 3 NSET SETRPZ-1', 
        odb=odb, 
        outputVariableName='Spatial displacement: U3 PI: rootAssembly Node 3 in NSET SETRPZ', 
        steps=('Step-1', ), __linkedVpName__='Viewport: 1')
    x0 = session.xyDataObjects['RF1 PI: rootAssembly N: 1 NSET SETRPX-1']
    x1 = session.xyDataObjects['RF1 PI: rootAssembly N: 2 NSET SETRPY-1']
    x2 = session.xyDataObjects['RF1 PI: rootAssembly N: 3 NSET SETRPZ-1']
    x3 = session.xyDataObjects['RF2 PI: rootAssembly N: 1 NSET SETRPX-1']
    x4 = session.xyDataObjects['RF2 PI: rootAssembly N: 2 NSET SETRPY-1']
    x5 = session.xyDataObjects['RF2 PI: rootAssembly N: 3 NSET SETRPZ-1']
    x6 = session.xyDataObjects['RF3 PI: rootAssembly N: 1 NSET SETRPX-1']
    x7 = session.xyDataObjects['RF3 PI: rootAssembly N: 2 NSET SETRPY-1']
    x8 = session.xyDataObjects['RF3 PI: rootAssembly N: 3 NSET SETRPZ-1']
    x9 = session.xyDataObjects['U1 PI: rootAssembly N: 1 NSET SETRPX-1']
    x10 = session.xyDataObjects['U1 PI: rootAssembly N: 2 NSET SETRPY-1']
    x11 = session.xyDataObjects['U1 PI: rootAssembly N: 3 NSET SETRPZ-1']
    x12 = session.xyDataObjects['U2 PI: rootAssembly N: 1 NSET SETRPX-1']
    x13 = session.xyDataObjects['U2 PI: rootAssembly N: 2 NSET SETRPY-1']
    x14 = session.xyDataObjects['U2 PI: rootAssembly N: 3 NSET SETRPZ-1']
    x15 = session.xyDataObjects['U3 PI: rootAssembly N: 1 NSET SETRPX-1']
    x16 = session.xyDataObjects['U3 PI: rootAssembly N: 2 NSET SETRPY-1']
    x17 = session.xyDataObjects['U3 PI: rootAssembly N: 3 NSET SETRPZ-1']
    session.xyReportOptions.setValues(numDigits=9, interpolation=ON)
    session.writeXYReport(
        fileName='C:/Users/aero_waas_2018_1/Documents/GitHub/composite-multiscale/micro_ruc/prepare_RUC_model/postprocessing/RUC_multi_strain_original.rpt', 
        appendMode=OFF, xyData=(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, 
        x11, x12, x13, x14, x15, x16, x17))


