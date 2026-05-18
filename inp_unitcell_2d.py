


#abaqus cae noGUI=C:\Users\p2515497\Documents\ghk-creep\ghk-creep/inp_unitcell_2d.py


# -*- coding: mbcs -*-
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *

import os 
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np

def ABQ_model(dt_loc):
	
	M = mdb.models['Model-1']
	
	# =============================================================================
	# Part
	# =============================================================================
	M.ConstrainedSketch(name='__profile__', sheetSize=200.0)
	M.sketches['__profile__'].rectangle(point1=(0.0, 0.0), 
	    point2=(1.0, 1.0))
	
	M.Part(dimensionality=TWO_D_PLANAR, name='Part-1', type=
	    DEFORMABLE_BODY)
	
	M.parts['Part-1'].BaseShell(sketch=
	    M.sketches['__profile__'])
	
	
	# =============================================================================
	# Material properties
	# =============================================================================
	
	M.Material(name='DP_NortonBailey')
	
	M.materials['DP_NortonBailey'].UserMaterial(
	    mechanicalConstants=(77e3, 0.2, 1.19, 0.0, 0.6, 3.77e-16, 0.716, 
	    -0.593))
	
	M.materials['DP_NortonBailey'].Depvar(n=5)
	
	M.HomogeneousSolidSection(material='DP_NortonBailey', name=
	    'Section-1', thickness=None)
	
	
	M.parts['Part-1'].SectionAssignment(offset=0.0, 
	    offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
	    faces=M.parts['Part-1'].faces.getSequenceFromMask(
	    mask=('[#1 ]', ), )), sectionName='Section-1', thicknessAssignment=
	    FROM_SECTION)
	M.rootAssembly.regenerate()
	# =============================================================================
	# Instance
	# =============================================================================
	
	I = M.rootAssembly.Instance(dependent=ON, name='Part-1-1', 
	    part=M.parts['Part-1'])
	
	# =============================================================================
	# Step 
	# =============================================================================
	M.ViscoStep(timePeriod=20000.0, initialInc=dt_loc, maxInc=dt_loc, maxNumInc=10000, minInc=0.001, 
				adaptiveDampingRatio=0.05, cetol=0.01, continueDampingFactors=False,
				name='Step-1', nlgeom=ON, previous='Initial', 
	    stabilizationMagnitude=0.0002, stabilizationMethod=DISSIPATED_ENERGY_FRACTION )
	
	
	# =============================================================================
	# Sets
	# =============================================================================
	M.rootAssembly.Set(edges=
	    I.edges.getSequenceFromMask(
	    ('[#1 ]', ), ), name='Btm')
	
	M.rootAssembly.Set(edges=
	   I.edges.getSequenceFromMask(
	    ('[#4 ]', ), ), name='TOP')
	
	M.rootAssembly.Set(name='RBM', vertices=
	   I.vertices.getSequenceFromMask(
	    ('[#1 ]', ), ))
	
	
	# =============================================================================
	#  BC
	# =============================================================================
	M.DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
	    distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
	    'BC-1', region=M.rootAssembly.sets['Btm'], u1=UNSET, 
	    u2=0.0, ur3=UNSET)
	
	
	M.rootAssembly.Surface(name='Surf-Top', side1Edges=
	   I.edges.getSequenceFromMask(
	    ('[#4 ]', ), ))
	M.Pressure(amplitude=UNSET, createStepName='Step-1', 
	    distributionType=UNIFORM, field='', magnitude=5.0, name='Load-1', region=
	    M.rootAssembly.surfaces['Surf-Top'])
	
	
	M.DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
	    distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
	    'BC-2', region=M.rootAssembly.sets['RBM'], u1=0.0, u2=
	    0.0, ur3=0.0)
	
	
	
	# =============================================================================
	# Mesh 
	# =============================================================================
	M.parts['Part-1'].seedPart(deviationFactor=0.1, 
	    minSizeFactor=0.1, size=1.0)
	
	M.parts['Part-1'].generateMesh()
	
	M.rootAssembly.regenerate()
	
	M.parts['Part-1'].setElementType(elemTypes=(ElemType(
	    elemCode=CPE4, elemLibrary=STANDARD), ElemType(elemCode=CPE3, 
	    elemLibrary=STANDARD)), regions=(
	    M.parts['Part-1'].faces.getSequenceFromMask(('[#1 ]', 
	    ), ), ))
	
	
	M.rootAssembly.Set(elements=
	   I.elements.getSequenceFromMask(
	    mask=('[#1 ]', ), ), name='ALL_LMT')
	
	
# 	M.keywordBlock.replace(23, 
# 	    '\n*Initial Conditions, type=SOLUTION\nALL_LMT, 0., 1.e-15, 0., 0., 0.')
	
	# =============================================================================
	# Output request
	# =============================================================================
	M.fieldOutputRequests['F-Output-1'].setValues(variables=(
	    'S', 'E', 'PE', 'PEEQ', 'PEMAG', 'CE', 'CEEQ', 'CEMAG', 'LE', 'U', 'RF', 
	    'CF', 'CSTRESS', 'CDISP', 'SDV'))
	# =============================================================================
	# Job
	# =============================================================================
	jobname = 'Creep_'+str(dt)
	mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
	    explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
	    memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, 
	    multiprocessingMode=DEFAULT, name=jobname, 
	    nodalOutputPrecision=SINGLE, numCpus=1, numGPUs=0, numThreadsPerMpiProcess=
	    1, queue=None, resultsFormat=ODB, scratch='', type=ANALYSIS, 
	    userSubroutine='', waitHours=0, waitMinutes=0)
	
	mdb.jobs[jobname].setValues(numThreadsPerMpiProcess=1, userSubroutine= '/home/duminy/Documents/02_Creep//umat_dp_norton_bailey.f')
	
	mdb.jobs[jobname].writeInput(consistencyChecking=OFF)
	

def PostProcess(odb):
	odbNew = str(path+odb[:-4])+'.odb'
	odb_object = openOdb(odbNew)
	session.viewports['Viewport: 1'].setValues(displayedObject=odb_object)
	V = (('SDV1', INTEGRATION_POINT),('SDV2', INTEGRATION_POINT),('SDV3', INTEGRATION_POINT),('SDV4', INTEGRATION_POINT))
	xyd = session.xyDataListFromField(odb = odb_object, outPosition=NODAL, variable = V, nodeSets = (RBM,))
	sdv1 = [[y for x,y in XY] for XY in xyd if XY.name=='SDV1']
	sdv2 = [[y for x,y in XY] for XY in xyd if XY.name=='SDV2']
	sdv3 = [[y for x,y in XY] for XY in xyd if XY.name=='SDV3']
	sdv4 = [[y for x,y in XY] for XY in xyd if XY.name=='SDV4']
	t 	 = [[x for x,y in XY] for XY in xyd if XY.name=='SDV4']

	return t,sdv1, sdv2, sdv3,sdv4


def Plot_Pq(p,q,titel):
	
	cm = 1 / 2.54
	fig, ax = plt.subplots(figsize=(8*cm, 7*cm))

	ax.plot((0,2),(0,2*np.tan(np.radians(71.57))),'-k')
	ax.plot((0,0),(0,2*np.tan(np.radians(71.57))),'-k')
	ax.plot((0,2*np.tan(np.radians(71.57))),(0,0),'-k')

	ax.text(2,2, 'Multi axial compression')
	ax.text(0.1,6.5, 'Shear-compression')

	points = [[0,0], [2,2*np.tan(np.radians(71.57))], [0, 6]]
	triangle = Polygon(points, fc=(1, 0, 0, 0.5), ec='red', lw=0)
	ax.add_patch(triangle)

	points = [[0,0], [2,2*np.tan(np.radians(71.57))], [6,0]]
	triangle2 = Polygon(points, fc=(0, 0, 1, 0.5), ec='blue', lw=0)
	ax.add_patch(triangle2)
	points = [[6,6], [2,2*np.tan(np.radians(71.57))], [6,0]]
	triangle2 = Polygon(points, fc=(0, 0, 1, 0.5), ec='blue', lw=0)
	ax.add_patch(triangle2)

	plt.arrow(0.5, 6.4, 0, -2, length_includes_head=True,
			head_width=0.1, head_length=0.1, color = 'red')

	plt.plot(p,q, 'ob')
	ax.set_aspect('equal', adjustable='box')
	ax.set_xlabel("p (MPa)")
	ax.set_ylabel("q (MPa)")
	plt.title(titel)


	plt.gca().set_frame_on(False)

	plt.show()



path = 'G:/01_Forschung/01_Creep/00_UMAT/UMAT_CLAUDE/UnitCell_2D/time_increment_dispControl/'
try:
	os.mkdir(path)
except WindowsError:
	pass

os.chdir(path)
delta_t = [1,5,25,150,500,750]

for dt in delta_t[:3]:
	#ABQ_model(dt)

	odb = 'Creep_'+str(dt)+'.odb'
	t,sd1,sd2,sd3,sd4 = PostProcess(odb)

	with open('Result_'+odb[:-4]+'.dat','w') as fid:
		for tti, s1,s2,s3,s4 in zip(t,sd1,sd2,sd3,sd4):
			fid.wirte(str(tti)+'\t'+str(s1)+'\t'+str(s2)'\t'+str(s3)+'\t'+str(-s4)+'\n')
	
	Plot_Pq(-sd4,sd3)

	
