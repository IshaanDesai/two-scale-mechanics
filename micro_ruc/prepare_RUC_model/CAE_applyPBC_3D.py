#
# Python script to implement periodic boundary conditions
#
# Christian Heinrich - cheinric@umich.edu
#

#ver 2012-01-22

def applyPBC():
##if 1:

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
    
    #from sys import path
    #path.append('C:\Python24\Lib\site-packages') 
    #from numpy import *
    
    modelName='Model-1'
    partName='Part-1'

    print 'Start \'clean up\''
    
    #'clean up'
    setNames=[]
    for key in mdb.models[modelName].parts[partName].sets.keys():
        if 'nodes' in key:
            setNames.append(key)
    setNames=tuple(setNames)          
    mdb.models['Model-1'].parts['Part-1'].deleteSets(setNames=setNames)

    constraintNames=[]          
    for key in mdb.models[modelName].constraints.keys():
        if 'PBCConstraint' in key:
            constraintNames.append(key)
    constraintNames=tuple(constraintNames)
    mdb.models['Model-1'].constraints.delete(constraintNames)
    print 'Done with \'clean up\''    

    
    #Tolerance for matching nodes
    tol=1e-5
   
    nodes={
        'x0':mdb.models[modelName].parts[partName].sets['setx0'].nodes,
        'x1':mdb.models[modelName].parts[partName].sets['setx1'].nodes,
        'y0':mdb.models[modelName].parts[partName].sets['sety0'].nodes,
        'y1':mdb.models[modelName].parts[partName].sets['sety1'].nodes,
        'z0':mdb.models[modelName].parts[partName].sets['setz0'].nodes,
        'z1':mdb.models[modelName].parts[partName].sets['setz1'].nodes}
    
    if not len(nodes['x0'])==len(nodes['x1']) or not len(nodes['y0'])==len(nodes['y1']) or not len(nodes['z0'])==len(nodes['z1']):
        print 'Warning: number of nodes is different'
        raise KeyboardInterrupt

    #check that bottom left is (0,0)
    allNodes=mdb.models[modelName].parts[partName].nodes
    xmin=1e9
    ymin=1e9
    zmin=1e9
    for i in range(len(allNodes)):
        coord=allNodes[i].coordinates
        xmin=min(xmin,coord[0])
        ymin=min(ymin,coord[1])
        zmin=min(zmin,coord[2])
    if xmin<-tol or ymin<-tol or zmin<-tol:
        print 'Warning: Corner not at (0,0,0)'
        raise KeyboardInterrupt        


    #Side length of the brick/ rectangle
    sideLength={'x':0.0,'y':0.0,'z':0.0}
    #For readibility...
    for node in nodes['x1']:
            sideLength['x']=max(node.coordinates[0],sideLength['x'])
    for node in nodes['y1']:
            sideLength['y']=max(node.coordinates[1],sideLength['y'])
    for node in nodes['z1']:
            sideLength['z']=max(node.coordinates[2],sideLength['z'])

    print 'lx='+str(sideLength['x'])+', ly='+str(sideLength['y'])+ ', lz='+ str(sideLength['z'])
    
    p=mdb.models[modelName].parts[partName]
    
    nodesLib={'x0':{},'x1':{},'y0':{},'y1':{},'z0':{},'z1':{}}
    
    #create the sets and note down their names and coordinates
    for key in nodes.keys():
        for i in range(len(nodes[key])):
            myNode=nodes[key][i:i+1]
            myLabel=nodes[key][i].label
            p.Set(nodes=myNode, name=key+'nodes'+str(myLabel))
            nodesLib[key][myLabel]=nodes[key][i].coordinates

    import sets
    # | union, & intersection, - difference
    #find the faces
    noEdge={}
    noEdge['x']=list((set(nodesLib['x0'].keys())|set(nodesLib['x1'].keys()))
        -(set(nodesLib['y0'].keys())|set(nodesLib['y1'].keys())|set(nodesLib['z0'].keys())|set(nodesLib['z1'].keys())))
    noEdge['y']=list((set(nodesLib['y0'].keys())|set(nodesLib['y1'].keys()))
        -(set(nodesLib['x0'].keys())|set(nodesLib['x1'].keys())|set(nodesLib['z0'].keys())|set(nodesLib['z1'].keys())))
    noEdge['z']=list((set(nodesLib['z0'].keys())|set(nodesLib['z1'].keys()))
        -(set(nodesLib['x0'].keys())|set(nodesLib['x1'].keys())|set(nodesLib['y0'].keys())|set(nodesLib['y1'].keys())))

    
    #find the corners
    corners={
        'x0y0z0':(),
        'x0y0z1':(),
        'x0y1z0':(),
        'x0y1z1':(),
        'x1y0z0':(),
        'x1y0z1':(),
        'x1y1z0':(),
        'x1y1z1':()}

    for key in corners.keys():
        corners[key]=list(set(nodesLib[key[0:2]].keys())&set(nodesLib[key[2:4]].keys())&set(nodesLib[key[4:6]].keys()))[0]

    #find the edges
    edges={
        'x0y0':(),
        'x0y1':(),
        'x1y0':(),
        'x1y1':(),
        'x0z0':(),
        'x0z1':(),
        'x1z0':(),
        'x1z1':(),
        'y0z0':(),
        'y0z1':(),
        'y1z0':(),
        'y1z1':()}

    for key in edges.keys():
        edges[key]=list((set(nodesLib[key[0:2]].keys())&set(nodesLib[key[2:4]].keys()))-set(corners.values()))
              
    pairs={'x':[],'y':[],'z':[]}

    #match the nodes
    for key in pairs.keys():
        for i in nodesLib[key+'0'].keys():
            if i in noEdge[key]:
                nodesLib1keys=nodesLib[key+'1'].keys()
                for j in nodesLib1keys:
                    diff=         abs(nodesLib[key+'1'][j][0]-nodesLib[key+'0'][i][0])
                    diff=diff+    abs(nodesLib[key+'1'][j][1]-nodesLib[key+'0'][i][1])
                    diff=diff+    abs(nodesLib[key+'1'][j][2]-nodesLib[key+'0'][i][2])
                    diff=abs(diff-sideLength[key])
                    if diff<tol:
                        pairs[key].append((i,j))
                        break
                    #When this is true, the code failed to find a partnering node
                    if j==nodesLib1keys[-1]:
                        print 'Could not find partnering node'
                        raise KeyboardInterrupt               

    
    constraintNumber=0        
            
    #create the equations on faces
    myKeys=pairs.keys()
    myKeys.sort()
    for key in myKeys:
        for i in range(len(pairs[key])):
            name1=partName+'-1.'+key+'1nodes'+str(pairs[key][i][1])
            name0=partName+'-1.'+key+'0nodes'+str(pairs[key][i][0])
            for j in [1,2,3]:
                mdb.models[modelName].Equation(name='PBCConstraint'+str(constraintNumber), terms=(
                    (1.0 ,name1, j), 
                    (-1.0,name0, j),
                    (-1.0,'setrp'+key,j)))
                constraintNumber=constraintNumber+1

    #create equations on corners
    #all nodes are either in x0 or xl
    for key in corners.keys():
        if corners[key] in nodesLib['x0'].keys():
            name=partName+'-1.x0nodes'+str(corners[key])            
        elif corners[key] in nodesLib['x1'].keys():
            name=partName+'-1.x1nodes'+str(corners[key])
        else:
            print 'Corner keys wrong!'
            raise KeyboardInterrupt
        
        for i in [1,2,3]:
            mdb.models[modelName].Equation(name='PBCConstraint'+str(constraintNumber), terms=(
                (1.0 ,name, i), 
                (-float(key[1]),'setrpx',i),
                (-float(key[3]),'setrpy',i),
                (-float(key[5]),'setrpz',i)))
            constraintNumber=constraintNumber+1

    #create equations for all edges
    for key in edges.keys():
        # If I am on a master edge, go to the beginning of the for-loop
        if key[1]=='0' and key[3]=='0':
            continue
        for nodeNumber in edges[key]:
            name1=partName+'-1.'+key[0:2]+'nodes'+str(nodeNumber)
            #find the corresponding node number in the zero side
            zeroEdgeNodes=edges[key[0]+'0'+key[2]+'0']
            if key[1]=='0':
                myNodesLib=nodesLib[key[2]+'1']
            else:
                myNodesLib=nodesLib[key[0]+'1']
                
            for zeroEdgeNode in zeroEdgeNodes:
                diff=abs(abs(myNodesLib[nodeNumber][0]-nodesLib[key[0]+'0'][zeroEdgeNode][0])
                        +abs(myNodesLib[nodeNumber][1]-nodesLib[key[0]+'0'][zeroEdgeNode][1])
                        +abs(myNodesLib[nodeNumber][2]-nodesLib[key[0]+'0'][zeroEdgeNode][2])
                        -sideLength[key[0]]*float(key[1])-sideLength[key[2]]*float(key[3]))
                if diff<tol:
                    name0=partName+'-1.'+key[0]+'0nodes'+str(zeroEdgeNode)
                    break
                if zeroEdgeNode==zeroEdgeNodes[-1]:
                    print 'Could not match edges nodes'
                    raise KeyboardInterrupt                    
            rpname0='setrp'+key[0]
            rpname1='setrp'+key[2]
            coeff0=-float(key[1])
            coeff1=-float(key[3])

            for j in [1,2,3]:
                mdb.models[modelName].Equation(name='PBCConstraint'+str(constraintNumber), terms=(
                    (1.0 ,name1, j),
                    (-1.0 ,name0, j),
                    (coeff0,rpname0,j),
                    (coeff1,rpname1,j)))
                constraintNumber=constraintNumber+1

    print '\n\n\nfinished applying periodic boundary conditions'
 
if __name__ == "__main__":
    applyPBC()  
    
    
