from abaqus import *
from abaqusConstants import *
import __main__

import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior







# units: um - s - kg - uN - MPa



# ##################################################################################
# geometric input values
# ##################################################################################


fiberRadius = (5.2/2.)/1000 # <<<
fiberVolumeFraction = 0.59  # <<<



lx = fiberRadius #this is through the thickness and somewhat arbitrary
towName = 'axial' # usually 'axial' or 'bias'

# other parameters
meshSize = 0.25/1000. 

deviationFactor = 0.05

# #towName='bias' # usually 'axial' or 'bias'
# #fiberVolumeFraction=0.77


# different loading cases and according boundary conditions to determine
# anisotropic material properties


loadingCases={}
loadingCases['s11']={'x1':0.1,'y3':0.0}
loadingCases['s22']={'y2':0.1,'x3':0.0}
loadingCases['s12']={'x2':0.1,'y1':0.0,'z2':0.0,'x3':0.0}
loadingCases['s23']={'y3':0.1,'z2':0.0,'x3':0.0,'y1':0.0}


# some intermidiate calculations

ly = sqrt(2*pi*fiberRadius**2./(sqrt(3)*fiberVolumeFraction))
lz = ly*sqrt(3)

lengths={'x':lx,'y':ly,'z':lz}

filename = 'out_axialSize.txt'
with open(filename, 'w') as f:
    f.write(str(lx) + ',' + str(ly) + ',' + str(lz) + '\n')
# close()


# Create new model database

Mdb()


# Create part


#extrusion
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=2*max(ly,lz))
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
s.rectangle(point1=(0.0, 0.0), point2=(lx,ly))
p = mdb.models['Model-1'].Part(name='Part-1', dimensionality=THREE_D, 
    type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Part-1']
p.BaseSolidExtrude(sketch=s, depth=lz)
s.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['Model-1'].sketches['__profile__']

#create sets for the side faces
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
faces = f.findAt(((0.0,ly/2.,lz/2.),),)
p.Set(faces=faces, name='setx0')
faces = f.findAt(((lx,ly/2.,lz/2.),),)
p.Set(faces=faces, name='setx1')
faces = f.findAt(((lx/2.,0.0,lz/2.),),)
p.Set(faces=faces, name='sety0')
faces = f.findAt(((lx/2.,ly,lz/2.),),)
p.Set(faces=faces, name='sety1')
faces = f.findAt(((lx/2.,ly/2.,0.0),),)
p.Set(faces=faces, name='setz0')
faces = f.findAt(((lx/2.,ly/2.,lz),),)
p.Set(faces=faces, name='setz1')

#partition face
p = mdb.models['Model-1'].parts['Part-1']
f, e, d1 = p.faces, p.edges, p.datums
t = p.MakeSketchTransform(sketchPlane=f[0], sketchUpEdge=e[0], 
    sketchPlaneSide=SIDE1, origin=(0.0, 0.0, 0.0))
s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
    sheetSize=2*max(ly,lz), gridSpacing=0.1, transform=t)
g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
s1.setPrimaryObject(option=SUPERIMPOSE)
p = mdb.models['Model-1'].parts['Part-1']
p.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)
s1.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(-fiberRadius, 0.0))
s1.CircleByCenterPerimeter(center=(0.0, ly), point1=(-fiberRadius, ly))
s1.CircleByCenterPerimeter(center=(lz, 0.0), point1=(lz+fiberRadius, 0.0))
s1.CircleByCenterPerimeter(center=(lz, ly), point1=(lz-fiberRadius, ly))
s1.CircleByCenterPerimeter(center=(lz/2., ly/2.), point1=(lz/2.+fiberRadius, ly/2.))
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedFaces = f.getSequenceFromMask(mask=('[#1 ]', ), )
e1, d2 = p.edges, p.datums
p.PartitionFaceBySketch(sketchUpEdge=e1[0], faces=pickedFaces, sketch=s1)
s1.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__profile__']

#extrude partitions
p = mdb.models['Model-1'].parts['Part-1']
c = p.cells
pickedCells = c[:]
e = p.edges
sweepPath =e.findAt(((lx/2.,0.0,0.0),),)[0]
pickedEdges=e.findAt(((0.0,fiberRadius*(sqrt(2.)/2.),fiberRadius*(sqrt(2.)/2.)),),)[0]
p.PartitionCellBySweepEdge(sweepPath=sweepPath, cells=pickedCells, 
        edges=(pickedEdges,))

pickedCells = c[:]
e = p.edges
sweepPath =e.findAt(((lx/2.,0.0,0.0),),)[0]
pickedEdges=e.findAt(((0.0,-fiberRadius*(sqrt(2.)/2.)+ly,fiberRadius*(sqrt(2.)/2.)),),)[0]
p.PartitionCellBySweepEdge(sweepPath=sweepPath, cells=pickedCells, 
        edges=(pickedEdges,))

pickedCells = c[:]
e = p.edges
sweepPath =e.findAt(((lx/2.,0.0,0.0),),)[0]
pickedEdges=e.findAt(((0.0,-fiberRadius*(sqrt(2.)/2.)+ly,-fiberRadius*(sqrt(2.)/2.)+lz),),)[0]
p.PartitionCellBySweepEdge(sweepPath=sweepPath, cells=pickedCells, 
        edges=(pickedEdges,))

pickedCells = c[:]
e = p.edges
sweepPath =e.findAt(((lx/2.,0.0,0.0),),)[0]
pickedEdges=e.findAt(((0.0,fiberRadius*(sqrt(2.)/2.),-fiberRadius*(sqrt(2.)/2.)+lz),),)[0]
p.PartitionCellBySweepEdge(sweepPath=sweepPath, cells=pickedCells, 
        edges=(pickedEdges,))

pickedCells = c[:]
e = p.edges
sweepPath =e.findAt(((lx/2.,0.0,0.0),),)[0]
pickedEdges=e.findAt(((0.0,-fiberRadius*(sqrt(2.)/2.)+ly/2.,-fiberRadius*(sqrt(2.)/2.)+lz/2.),),)[0]
p.PartitionCellBySweepEdge(sweepPath=sweepPath, cells=pickedCells, 
        edges=(pickedEdges,))

#create sets
p = mdb.models['Model-1'].parts['Part-1']
c = p.cells
cells = c.findAt(((lx/2.,0.,0.),),((lx/2.,ly,lz),),((lx/2.,ly,0.),),((lx/2.,0.,lz),),((lx/2.,ly/2.,lz/2.),))
p.Set(cells=cells, name='setFiber')
cells = c.findAt(((lx/2.,ly/2.,0.),),)
p.Set(cells=cells, name='setMatrix')

#create coordinate system
p.DatumCsysByThreePoints(name='Datum csys-1', coordSysType=CARTESIAN, origin=(
    0.0, 0.0, 0.0), line1=(1.0, 0.0, 0.0), line2=(0.0, 1.0, 0.0))


# Mesh part


p = mdb.models['Model-1'].parts['Part-1']
c = p.cells

#some more partitions for better meshing

pickedCells = p.cells[:]
e1, v1, d2 = p.edges, p.vertices, p.datums
p.PartitionCellByPlanePointNormal(normal=e1.findAt(((0,0,lz/2.),),)[0], cells=pickedCells, 
        point=p.InterestingPoint(edge=e1.findAt(((0,0,lz/2.),),)[0], rule=MIDDLE))

pickedCells = p.cells[:]
e1, v1, d2 = p.edges, p.vertices, p.datums
p.PartitionCellByPlanePointNormal(normal=e1.findAt(((0,ly/2.,0),),)[0], cells=pickedCells, 
        point=p.InterestingPoint(edge=e1.findAt(((0,ly/2.,0),),)[0], rule=MIDDLE))

#assign meshing technique: swept mesh
c = p.cells
pickedRegions = c.getSequenceFromMask(mask=('[#fff ]', ), )
p.setMeshControls(regions=pickedRegions, technique=SWEEP, 
    algorithm=ADVANCING_FRONT)

#set element type
elemType1 = mesh.ElemType(elemCode=C3D20, elemLibrary=STANDARD)
elemType2 = mesh.ElemType(elemCode=C3D15, elemLibrary=STANDARD)
elemType3 = mesh.ElemType(elemCode=C3D10, elemLibrary=STANDARD)
pickedRegions =(c[:], )
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, 
    elemType3))

#seed
p.seedPart(size=meshSize, deviationFactor=deviationFactor, minSizeFactor=0.1)

#seed fiber radii
e=mdb.models['Model-1'].parts['Part-1'].edges
for i in range(len(e)):
    if abs(e[i].getSize()-fiberRadius)<1e-6:
        print e[i].getSize()
        pickedEdges = e[i:i+1]
        p.seedEdgeBySize(edges=pickedEdges, size=meshSize*4., deviationFactor=0.1, 
            constraint=FINER)

#mesh cells piece by piece to ensure that node match up
pickedRegions=p.sets['setFiber'].cells[:]
p.generateMesh(regions=pickedRegions)

pickedRegions=p.sets['setMatrix'].cells[:]
p.generateMesh(regions=pickedRegions)


# Create materials


mdb.models['Model-1'].Material(name='materialfiber')
mdb.models['Model-1'].Material(name='materialMatrix')

#get the input values from file
#fiber
if os.path.isfile('input_fiberProps.txt'):
    f=open('input_fiberProps.txt')
    lines=f.readlines()
    f.close()
    line=lines[2].strip().split(',')
    engineeringConstFiber=[0,0,0,0,0,0,0,0,0]
    engineeringConstFiber[0]=float(line[0])#/1e6
    engineeringConstFiber[1]=float(line[1])#/1e6
    engineeringConstFiber[2]=float(line[2])#/1e6
    engineeringConstFiber[3]=float(line[3])
    engineeringConstFiber[4]=float(line[4])
    engineeringConstFiber[5]=float(line[5])
    engineeringConstFiber[6]=float(line[6])#/1e6
    engineeringConstFiber[7]=float(line[7])#/1e6
    engineeringConstFiber[8]=float(line[8])#/1e6
    engineeringConstFiber=tuple(engineeringConstFiber)
else:
    print 'Could not find files with fiber properties'
    raise

#matrix
#import the material properties from file
#axial properties
if os.path.isfile('input_matrixProps.txt'):
    f=open('input_matrixProps.txt')
    lines=f.readlines()
    f.close()
    line=lines[2].strip().split(',')
    engineeringConstMatrix=[0,0]
    engineeringConstMatrix[0]=float(line[0])#/1e6
    engineeringConstMatrix[1]=float(line[1])
    engineeringConstMatrix=tuple(engineeringConstMatrix)
    yieldTableMatrix=[]
    for i in range(len(lines)-6):
        yieldTableMatrix.append((float(lines[i+6].strip().split(',')[0]),float(lines[i+6].strip().split(',')[1])))
        # yieldTableMatrix.append((float(lines[i+6].strip().split(',')[0])/1e6,float(lines[i+6].strip().split(',')[1])))
    yieldTableMatrix=tuple(yieldTableMatrix)
else:
    print 'Could not find files with matrix properties'
    raise

#assign material values to material properties
mdb.models['Model-1'].materials['materialfiber'].Elastic(type=ENGINEERING_CONSTANTS, table=(engineeringConstFiber, ))
mdb.models['Model-1'].materials['materialMatrix'].Elastic(table=(engineeringConstMatrix, ))
mdb.models['Model-1'].materials['materialMatrix'].Plastic(table=yieldTableMatrix)

mdb.models['Model-1'].HomogeneousSolidSection(name='sectionFiber', 
    material='materialfiber', thickness=None)
mdb.models['Model-1'].HomogeneousSolidSection(name='sectionMatrix', 
    material='materialMatrix', thickness=None)


# Create and assign sections


mdb.models['Model-1'].HomogeneousSolidSection(name='sectionFiber', 
    material='materialFiber', thickness=None)
mdb.models['Model-1'].HomogeneousSolidSection(name='sectionMatrix', 
    material='materialMatrix', thickness=None)

p = mdb.models['Model-1'].parts['Part-1']

region = p.sets['setFiber']
p.SectionAssignment(region=region, sectionName='sectionFiber', offset=0.0, 
    offsetType=MIDDLE_SURFACE, offsetField='', 
    thicknessAssignment=FROM_SECTION)

region = p.sets['setMatrix']
p.SectionAssignment(region=region, sectionName='sectionMatrix', offset=0.0, 
    offsetType=MIDDLE_SURFACE, offsetField='', 
    thicknessAssignment=FROM_SECTION)

#assign material orientation
p = mdb.models['Model-1'].parts['Part-1']
region = p.sets['setFiber']
orientation = mdb.models['Model-1'].parts['Part-1'].datums.items()[0][1]
mdb.models['Model-1'].parts['Part-1'].MaterialOrientation(region=region, 
    orientationType=SYSTEM, axis=AXIS_3, localCsys=orientation, 
    fieldName='', additionalRotationType=ROTATION_NONE, angle=0.0, 
    additionalRotationField='', stackDirection=STACK_3)


# Create assembly

#copy part
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = mdb.models['Model-1'].parts['Part-1']
a.Instance(name='Part-1-1', part=p, dependent=ON)

#create reference points
a.ReferencePoint(point=(lx, 0.0, 0.0))
a.ReferencePoint(point=(0.0, ly, 0.0))
a.ReferencePoint(point=(0.0, 0.0, lz))

#create reference point sets
a.Set(referencePoints=((a.referencePoints.findAt((lx,0.0,0.0),),)), name='setrpx')
a.Set(referencePoints=((a.referencePoints.findAt((0.0,ly,0.0),),)), name='setrpy')
a.Set(referencePoints=((a.referencePoints.findAt((0.0,0.0,lz),),)), name='setrpz')

# Create step
mdb.models['Model-1'].StaticStep(name='Step-1', previous='Initial', 
    timePeriod=1., maxNumInc=10000, initialInc=1e-7, minInc=1e-7, 
    maxInc=0.1, nlgeom=ON)

# Create field output
del mdb.models['Model-1'].fieldOutputRequests['F-Output-1']

#mdb.models['Model-1'].TimePoint(name='TimePoints-1', points=((0.0, 0.5, 0.05), (0.5, 1.0, 0.1)))
mdb.models['Model-1'].TimePoint(name='TimePoints-1', points=((0.5, ), (1.0, )))

mdb.models['Model-1'].FieldOutputRequest(name='F-Output-1', 
    createStepName='Step-1', variables=('S', 'MISESMAX', 'E', 'U'), 
    timePoint='TimePoints-1')




# History output


del mdb.models['Model-1'].historyOutputRequests['H-Output-1']

regionDef=mdb.models['Model-1'].rootAssembly.sets['setrpx']
mdb.models['Model-1'].HistoryOutputRequest(name='H-Output-1', 
    createStepName='Step-1', variables=('U1', 'U2', 'U3', 'RF1', 'RF2', 
    'RF3'), region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)
mdb.models['Model-1'].historyOutputRequests['H-Output-1'].setValues(
    numIntervals=500)

regionDef=mdb.models['Model-1'].rootAssembly.sets['setrpy']
mdb.models['Model-1'].HistoryOutputRequest(name='H-Output-2', 
    createStepName='Step-1', variables=('U1', 'U2', 'U3', 'RF1', 'RF2', 
    'RF3'), region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)
mdb.models['Model-1'].historyOutputRequests['H-Output-2'].setValues(
    numIntervals=500)

regionDef=mdb.models['Model-1'].rootAssembly.sets['setrpz']
mdb.models['Model-1'].HistoryOutputRequest(name='H-Output-3', 
    createStepName='Step-1', variables=('U1', 'U2', 'U3', 'RF1', 'RF2', 
    'RF3'), region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)
mdb.models['Model-1'].historyOutputRequests['H-Output-3'].setValues(
    numIntervals=500)


# Apply periodic boundary conditions


import CAE_applyPBC_3D
CAE_applyPBC_3D.applyPBC()


# consider the 4 different cases to find transversly isotropic peoperties


for loadingCase in loadingCases:
    #make a new model for each
    modelName='Model-'+loadingCase
    mdb.Model(name=modelName, objectToCopy=mdb.models['Model-1'])
    

    # Assign boundary conditions

    a = mdb.models[modelName].rootAssembly

    #delete all previous boundary conditions
    for item in mdb.models[modelName].boundaryConditions.items():
       del mdb.models[modelName].boundaryConditions[item[0]]

    for key in loadingCases[loadingCase]:
        u={'1':UNSET,'2':UNSET,'3':UNSET}
        setrp='setrp'+key[0]
        u[key[1]]=loadingCases[loadingCase][key]*lengths[key[0]]

        region = a.sets[setrp]
        mdb.models[modelName].DisplacementBC(name='BC-'+key, createStepName='Step-1', 
            region=region, u1=u['1'], u2=u['2'], u3=u['3'], ur1=UNSET, ur2=UNSET, 
            ur3=UNSET, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, 
            fieldName='', localCsys=None)

    # ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
    # Create job and input file
    # ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    jobName='Job_'+towName+'_'+loadingCase

    mdb.Job(name=jobName, model=modelName, description='', type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=75, 
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
        scratch='', multiprocessingMode=DEFAULT, numCpus=4, numDomains=4)

    mdb.jobs[jobName].writeInput(consistencyChecking=OFF)

    
