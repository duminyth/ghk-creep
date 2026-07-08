


#abaqus cae noGUI=C:\Users\p2515497\Documents\ghk-creep\ghk-creep\inp_unitcell_2d.py


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
	V = (('SDV1', INTEGRATION_POINT),('SDV2', INTEGRATION_POINT),('SDV3', INTEGRATION_POINT),('SDV4', INTEGRATION_POINT),('SDV5', INTEGRATION_POINT), ('LE', INTEGRATION_POINT, ((COMPONENT, 'LE22'),)))
	
	#xyd = session.xyDataListFromField(odb = odb_object, outputPosition=ELEMENT_NODAL, variable = V, nodeSets = ('RBM',))
	xyd = session.xyDataListFromField(odb = odb_object, outputPosition=ELEMENT_NODAL, variable = V, nodeSets = ('PART-1-1.RBMC',))
    
	sdv1 = [[y for x, y in XY] for XY in xyd if 'SDV1' in XY.name]
	sdv2 = [[y for x, y in XY] for XY in xyd if 'SDV2' in XY.name]
	sdv3 = [[y for x, y in XY] for XY in xyd if 'SDV3' in XY.name]
	sdv4 = [[y for x, y in XY] for XY in xyd if 'SDV4' in XY.name]
	sdv5 = [[y for x, y in XY] for XY in xyd if 'SDV5' in XY.name]
	LE 	 = [[y for x, y in XY] for XY in xyd if 'LE22' in   XY.name]
	t    = [[x for x, y in XY] for XY in xyd if 'SDV4' in XY.name]
    
	return t,sdv1, sdv2, sdv3,sdv4, sdv5, LE


def Plot_Pq(p,q,titel):
    
	cm = 1 / 2.54
	fig, ax = plt.subplots(figsize=(8*cm, 7*cm))

	p = np.asarray(p)
	q = np.asarray(q)
	maxX = np.max(p)
	maxq = np.max(q)


	angle = np.radians(71.57)
	p_yield_top = maxq / np.tan(angle)

	ax.plot((0,0),(0,maxq*1.1),'-k')
	ax.plot((0,0),(0,maxX*1.1),'-k')
	ax.plot((maxq/np.tan(np.radians(71.57)),maxq),(0,0),'-k')
	
	if True:
		points = [[0,0], [maxq/np.tan(np.radians(71.57)),maxq], [0, maxq]]
		triangle = Polygon(points, fc=(1, 0, 0, 0.5), ec='red', lw=0)
		ax.add_patch(triangle)
		
		points = [[0,0], [maxq/np.tan(np.radians(71.57)),maxq], [6,0]]
		triangle2 = Polygon(points, fc=(0, 0, 1, 0.5), ec='blue', lw=0)
		ax.add_patch(triangle2)
		points = [[maxX,maxq], [maxq/np.tan(np.radians(71.57)),maxq], [6,0]]
		triangle2 = Polygon(points, fc=(0, 0, 1, 0.5), ec='blue', lw=0)
		ax.add_patch(triangle2)
		
		
	plt.plot(p,q, 'ob')
	ax.set_aspect('equal', adjustable='box')
	ax.set_xlabel("p (MPa)")
	ax.set_ylabel("q (MPa)")
    
	ax.plot((0,0),(0,maxq*1.1),'-k')
	ax.plot((0,0),(0,maxX*1.1),'-k')
	ax.plot((maxq/np.tan(np.radians(71.57)),maxq),(0,0),'-k')
    
	plt.gca().set_frame_on(False)
	plt.savefig(path+str(titel)+'.png', dpi=300)
    
	
	plt.close(fig)

	cm = 1 / 2.54
	fig, ax = plt.subplots(figsize=(8*cm, 7*cm))
	maxX = np.max(p)
	maxq = np.max(q)
    
	ax.plot((0,0),(0,maxq*1.1),'-k')
	ax.plot((0,0),(0,maxX*1.1),'-k')
	ax.plot((maxq/np.tan(np.radians(71.57)),maxq),(0,0),'-k')
    
    
	points = [[0,0], [maxq/np.tan(np.radians(71.57)),maxq], [0, maxq]]
	triangle = Polygon(points, fc=(1, 0, 0, 0.5), ec='red', lw=0)
	ax.add_patch(triangle)
    
	points = [[0,0], [maxq/np.tan(np.radians(71.57)),maxq], [6,0]]
	triangle2 = Polygon(points, fc=(0, 0, 1, 0.5), ec='blue', lw=0)
	ax.add_patch(triangle2)
	points = [[maxX,maxq], [maxq/np.tan(np.radians(71.57)),maxq], [6,0]]
	triangle2 = Polygon(points, fc=(0, 0, 1, 0.5), ec='blue', lw=0)
	ax.add_patch(triangle2)
    
	plt.arrow(0.5, 1.15*maxq, 0, -2, length_includes_head=True,
			head_width=0.1, head_length=0.1, color = 'red')
    
	plt.plot(p,q, 'ob')
	ax.set_aspect('equal', adjustable='box')
	ax.set_aspect('auto')
	ax.set_xlabel("p (MPa)")
	ax.set_ylabel("q (MPa)")
    
	ax.plot((0,0), (0, maxq*1.1), '-k')
	ax.plot((0,maxX*1.1), (0,0), '-k')
	ax.plot((maxq/np.tan(np.radians(71.57)),maxq),(0,0),'-k')
	ax.set_ylim(5,15)
	ax.set_xlim(2.5,9)
    
	plt.gca().set_frame_on(False)
	plt.savefig(path+str(titel)+'_zoomed.png', dpi=300)
    
	
	plt.close(fig)
	return


from matplotlib.patches import Polygon

def Plot_Pq2(p, q, titel, show_zones=True):
    
    cm = 1 / 2.54
    fig, ax = plt.subplots(figsize=(10*cm, 7*cm))

    p = np.asarray(p)
    q = np.asarray(q)

    maxp = np.max(p)
    maxq = np.max(q)

    angle = np.radians(71.57)
    p_yield_top = maxq / np.tan(angle)

    # Colored zones
    if show_zones:
        points_red = [
            [0, 0],
            [p_yield_top, maxq],
            [0, maxq]
        ]
        triangle_red = Polygon(
            points_red,
            facecolor=(1, 0, 0, 0.25),
            edgecolor='none'
        )
        ax.add_patch(triangle_red)

        points_blue = [
            [0, 0],
            [p_yield_top, maxq],
            [maxp, 0]
        ]
        triangle_blue = Polygon(
            points_blue,
            facecolor=(0, 0, 1, 0.20),
            edgecolor='none'
        )
        ax.add_patch(triangle_blue)

    # Axes
    ax.plot([0, 0], [0, maxq*1.1], '-k', lw=1.0)
    ax.plot([0, maxp*1.1], [0, 0], '-k', lw=1.0)

    # Yield line
    ax.plot([0, p_yield_top], [0, maxq], '-k', lw=1.2)

    # Data
    ax.plot(p, q, '-ob', markersize=2.0, linewidth=1.0)

    ax.set_xlabel("p (MPa)")
    ax.set_ylabel("q (MPa)")

    ax.set_xlim(0, maxp*1.05)
    ax.set_ylim(0, maxq*1.15)

    ax.set_aspect('auto')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.tight_layout()
    fig.savefig(path + str(titel) + '.png', dpi=300)
    plt.close(fig)


def Plot_eps(e_pl,e_cr,titel):
    
	cm = 1 / 2.54
	fig, ax = plt.subplots(figsize=(8*cm, 7*cm))
    
        
	plt.plot(e_pl,e_cr, 'ob')
	ax.set_xlabel(r"$\varepsilon_\mathrm{pl}$")
	ax.set_ylabel(r"$\varepsilon_\mathrm{cr}$")
	plt.title(titel)
    
	ax.set_aspect('auto', adjustable='box')
	plt.gca().set_frame_on(True)
	plt.gca().ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
	plt.tight_layout()
    
    
    
	plt.savefig(path+'eps_dt'+str(dt)+'.png', dpi=300, bbox_inches="tight")
	
	plt.close(fig)
	return

def Plot_check_0(t,eps_pl, eps_cr, q, p, Deps_cr, titel):
	cm = 1 / 2.54
	fig, axes = plt.subplots(
		4, 1,
		sharex=True,
		figsize=(8*cm, 8*cm)
	)
    
	axes[0].plot(t, eps_pl)
	axes[0].set_ylabel(r"$\varepsilon_\mathrm{pl}$")
    
	axes[1].plot(t, eps_cr)
	axes[1].set_ylabel(r"$\varepsilon_\mathrm{cr}$")
    
	axes[2].plot(t, q)
	axes[2].set_ylabel("q (MPa)")
    
	axes[3].plot(t, -np.array(p))
	axes[3].set_ylabel("p(MPa)")
	axes[3].set_xlabel("time(s)")
    
	for ax in axes:
		ax.grid(True)
    
	plt.tight_layout()
	plt.savefig(path+'Evolution_qty'+str(dt)+'_'+str(titel)+'.png', dpi=300)
    
	
	plt.close(fig)
    
	cm = 1 / 2.54
	fig, ax = plt.subplots(figsize=(8*cm, 7*cm))
	ax.plot(t, Deps_cr)
	ax.set_ylabel(r"$\Delta\varepsilon_\mathrm{cr}$")
	ax.set_xlabel("$time(s)$")
    
	ax.set_aspect('auto', adjustable='box')
	plt.gca().set_frame_on(True)
	plt.tight_layout()
	plt.savefig(path+'DelatCreep'+str(dt)+'_'+str(titel)+'.png', dpi=300,bbox_inches="tight")
	
	plt.close(fig)
    
    
    
	cm = 1 / 2.54
	fig, ax = plt.subplots(figsize=(8*cm, 7*cm))
	ax.plot(t, eps_cr)
	ax.set_ylabel(r"$\varepsilon_\mathrm{cr}$")
	ax.set_xlabel("$time(s)$")
    
	ax.set_aspect('auto', adjustable='box')
	plt.gca().set_frame_on(True)
	plt.tight_layout()
	plt.gca().ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
	plt.savefig(path+'Eps_Creep'+str(dt)+'_'+str(titel)+'.png', dpi=300,bbox_inches="tight")
	
	plt.close(fig)
    
	return


def Plot_stress_all(results):
	cm = 1 / 2.54
	fig, axes = plt.subplots(
		2, 1,
		sharex=True,
		figsize=(11*cm, 7*cm)
	)
	for res in results:
		axes[0].plot(res[0], res[1], label=res[3])
		axes[0].set_ylabel("q (MPa)")
        
		axes[1].plot(res[0], -np.array(res[2]), label=res[3])
		axes[1].set_ylabel("p (MPa)")
		axes[1].set_xlabel("time (s)")
        
		axes[1].legend()
	plt.tight_layout()
	
	plt.savefig(path+'Evolution_pq_all.png', dpi=300)
	
	plt.close(fig)
    
    
	
    
	fig, axes = plt.subplots(figsize=(11*cm, 4*cm))
    
	for res in results:
		strain = np.array(res[4])+np.array(res[5])+np.array(res[6])
		axes.plot(res[0], strain, label=res[3])
		axes.set_ylabel(r"$\varepsilon_\mathrm{tot}$")
        
		axes.set_xlabel("time (s)")
        
		axes.legend()
	plt.tight_layout()
	
	plt.savefig(path+'strain_all.png', dpi=300)
	
	plt.close(fig)
    
	return

def Extract_pq_ABQ_noUMAT(odb, path):
	odbNew = str(path+odb[:-4])+'.odb'
	odb_object = openOdb(odbNew)
	session.viewports['Viewport: 1'].setValues(displayedObject=odb_object)
	V = (('S', INTEGRATION_POINT),)
	
	xyd = session.xyDataListFromField(odb = odb_object, outputPosition=ELEMENT_NODAL, variable = V, nodeSets = ('RBM',))
    

	VAR_S11   = np.array([[y for x,y in XY] for XY in xyd if 'S11' in XY.name])[0]
	VAR_S22   = np.array([[y for x,y in XY] for XY in xyd if 'S22' in XY.name])[0]
	VAR_S33   = np.array([[y for x,y in XY] for XY in xyd if 'S33' in XY.name])[0]
	VAR_S12   = np.array([[y for x,y in XY] for XY in xyd if 'S12' in XY.name])[0]
	t   = np.array([[x for x,y in XY] for XY in xyd if 'S12' in XY.name])[0]


	VAR_S13   = 0.
	VAR_S23   = 0.
    
	p = (VAR_S11+VAR_S22+VAR_S33)*1/3.
    
	q = np.sqrt(3/2)*np.sqrt(VAR_S11**2+VAR_S22**2+VAR_S12**2+VAR_S13**2+VAR_S23**2)
    
	return p,q,t


def Extract_Strain_ABQ_noUMAT(odb, path):
	odbNew = str(path+odb[:-4])+'.odb'
	odb_object = openOdb(odbNew)
	session.viewports['Viewport: 1'].setValues(displayedObject=odb_object)
	V = (('CEEQ', INTEGRATION_POINT),)
	
	xyd = session.xyDataListFromField(odb = odb_object, outputPosition=ELEMENT_NODAL, variable = V, nodeSets = ('RBM',))
    

	VAR_CEEQ   = np.array([[y for x,y in XY] for XY in xyd if 'CEEQ' in XY.name])[0]

    
	return VAR_CEEQ


def Extract_Strain_ABQ_CreepRoutine(odb, path):
	odbNew = str(path+odb[:-4])+'.odb'
	odb_object = openOdb(odbNew)
	session.viewports['Viewport: 1'].setValues(displayedObject=odb_object)
	V = (('SDV1', INTEGRATION_POINT),)
	
	xyd = session.xyDataListFromField(odb = odb_object, outputPosition=ELEMENT_NODAL, variable = V, nodeSets = ('RBM',))
    

	VAR_CEEQ   	= [[y for x, y in XY] for XY in xyd if 'SDV1' in XY.name]
	VAR_CEEQ   	= np.array([[y for x,y in XY] for XY in xyd if 'SDV1' in XY.name])[0]
	t   		= [[x for x, y in XY] for XY in xyd if 'SDV1' in XY.name][0]

	cm = 1 / 2.54
	fig, ax = plt.subplots(figsize=(8*cm, 4*cm))
	ax.plot(t, VAR_CEEQ)
	ax.set_ylabel(r"$\varepsilon_\mathrm{cr}$")
	ax.set_xlabel("$time(s)$")
    
	ax.set_aspect('auto', adjustable='box')
	plt.gca().set_frame_on(True)
	plt.tight_layout()
	plt.gca().ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
	plt.savefig(path+'Eps_Creep'+str(dt)+'.png', dpi=300,bbox_inches="tight")
	
	plt.close(fig)
    
	return VAR_CEEQ


def Extract_Strain_ABQ_Nothing(odb, path, name):
	odbNew = str(path+odb[:-4])+'.odb'
	odb_object = openOdb(odbNew)
	session.viewports['Viewport: 1'].setValues(displayedObject=odb_object)
	V = (('CEEQ', INTEGRATION_POINT),('PEEQ', INTEGRATION_POINT),('EE', INTEGRATION_POINT),('CE', INTEGRATION_POINT),('PE', INTEGRATION_POINT),)
	
	xyd = session.xyDataListFromField(odb = odb_object, outputPosition=ELEMENT_NODAL, variable = V, nodeSets = ('RBM',))
    

	VAR_CEEQ   	= np.array([[y for x,y in XY] for XY in xyd if 'CEEQ' in XY.name])[0]
	VAR_PEEQ   	= np.array([[y for x,y in XY] for XY in xyd if 'PEEQ' in XY.name])[0]
	VAR_EE22   	= np.array([[y for x,y in XY] for XY in xyd if 'EE22' in XY.name])[0]
	VAR_CE22   	= np.array([[y for x,y in XY] for XY in xyd if 'CE22' in XY.name])[0]
	VAR_PE22   	= np.array([[y for x,y in XY] for XY in xyd if 'PE22' in XY.name])[0]
	t   		= [[x for x, y in XY] for XY in xyd if 'CEEQ' in XY.name][0]

	cm = 1 / 2.54
	fig, ax = plt.subplots(figsize=(8*cm, 4*cm))
	ax.plot(t, VAR_CEEQ)
	ax.set_ylabel(r"$\varepsilon_\mathrm{cr}$")
	ax.set_xlabel("$time(s)$")
    
	ax.set_aspect('auto', adjustable='box')
	plt.gca().set_frame_on(True)
	plt.tight_layout()
	plt.gca().ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
	plt.savefig(path+'Eps_Creep'+str(dt)+'_'+str(name)+'.png', dpi=300,bbox_inches="tight")
	
	plt.close(fig)

	cm = 1 / 2.54
	fig, ax = plt.subplots(figsize=(8*cm, 4*cm))
	ax.plot(t, VAR_PEEQ)
	ax.set_ylabel(r"$\varepsilon_\mathrm{pl}$")
	ax.set_xlabel("$time(s)$")
    
	ax.set_aspect('auto', adjustable='box')
	plt.gca().set_frame_on(True)
	plt.tight_layout()
	plt.gca().ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
	plt.savefig(path+'Eps_Plas'+str(dt)+'_'+str(name)+'.png', dpi=300,bbox_inches="tight")
	
	plt.close(fig)
    
	return VAR_CEEQ,VAR_PEEQ, VAR_EE22, VAR_CE22, VAR_PE22


def Plot_stress_strain(p,q,ce,pe,ee, titel):
    
	cm = 1 / 2.54
	fig, ax = plt.subplots(figsize=(8*cm, 7*cm))
	maxX = np.max(p)
	maxq = np.max(q)
    
	total_strain = np.abs(ce+pe+ee)
    
	plt.plot(total_strain,q, 'ob', label='q')
	plt.plot(total_strain,p, '+r', label='p')
	ax.set_aspect('equal', adjustable='box')
	ax.set_xlabel(r"$\varepsilon_\mathrm{pl}$")
	ax.set_ylabel(r"$\sigma$")
	ax.legend()

	ax.set_aspect('auto', adjustable='box')
	plt.gca().set_frame_on(True)
	plt.tight_layout()
	plt.savefig(path+str(titel)+'.png', dpi=300)
    
	
	plt.close(fig)

	with open(path+'stress_strain_data_'+str(titel)+'.dat', 'w') as fid:
		for x,y,z in zip(total_strain,p,q):
			fid.write(str(x)+'\t'+str(y)+'\t'+str(z)+'\n')
	return


path = 'G:/01_Forschung/01_Creep/00_UMAT/UMAT_CLAUDE/UnitCell_2D/time_increment_umat_v2/1E-13/'
path = 'G:/01_Forschung/01_Creep/00_UMAT/UMAT_CLAUDE/UnitCell_2D/time_increment_umat_v2/Effect_decr_sub/K=1000/'
#path = 'G:/01_Forschung/01_Creep/00_UMAT/UMAT_CLAUDE/UnitCell_2D/time_increment_umat_v2/Effect_ecr0/K_1000/'
path = 'G:/01_Forschung/01_Creep/00_UMAT/UMAT_CLAUDE/Comparing_Method/NewUMAT/CompressionShear/'
path = 'G:/01_Forschung/01_Creep/00_UMAT/UMAT_CLAUDE/UnitCell3D/'

try:
	os.mkdir(path)
except WindowsError:
	pass

os.chdir(path)
delta_t = [1,5,25,150,500,750]
delta_t = [2,6,10,12,15]
delta_t = [4,10,12,14]
p_all = []
q_all = []
if True:
	for dt in delta_t[1:2]:

		if True:
			#New UMAT
				#ABQ_model(dt)
				
				odb = 'E-12.odb'
				odb = 'k_2_5_1p_cs.odb'
				odb = 'UniCell3D_Clamped_1C3S.odb'
				odb = 'TwentySevenCell_3D.odb'
				#odb = 'Creep_'+str(dt)+'.odb'
				t,sd1,sd2,sd3,sd4,sd5,LE = PostProcess(odb)
				
				with open('Result_'+odb[:-4]+'.dat','w') as fid:
					for tti, s1,s2,s3,s4,s5 in zip(t[0],sd1[0],sd2[0],sd3[0],sd4[0], sd5[0]):
						fid.write(str(tti)+"\t"+str(s1)+"\t"+str(s2)+'\t'+str(s3)+'\t'+str(-s4)+'\t'+str(s5)+'\n')
				
				p_all.append((t[0], sd3[0],sd4[0], dt,LE[0],sd1[0],sd2[0]))
				
				
				Plot_Pq2(-np.array(sd4[0]),sd3[0], 'TwentySevenCell_3D')
				Plot_Pq2(-np.array(sd4[0][:65]),sd3[0][:65], 'TwentySevenCell_3D_zoom')
				Plot_eps(sd1[0],sd2[0], "")
				Plot_check_0(t[0],sd1[0],sd2[0],sd3[0],sd4[0], sd5[0], 'TwentySevenCell_3D')
	if False:
		Plot_stress_all(p_all)
    
	if True:
		odb = 'G1_schachter.odb'
		t,sd1,sd2,sd3,sd4,sd5,LE = PostProcess(odb)
    
	# ABQ model 

	if False:
		odb = 'Creep_Abq_0_5_noDamp_CompressionShear.odb'
		path = 'G:/01_Forschung/01_Creep/00_UMAT/UMAT_CLAUDE/Comparing_Method/UnitCell_Comparison_Jin/'
		p,q,t = Extract_pq_ABQ_noUMAT(odb, path )
		ceeq, peeq, ee,ce,pe = Extract_Strain_ABQ_Nothing(odb, path, '1_P_noDamp_CompressionShear')
		Plot_Pq2(-np.array(p),q, 'alpha_0_5_1P_noDamp_CompressionShear')
		Plot_Pq2(-np.array(p)[:65],q[:65], 'alpha_0_5_1P_noDamp_CompressionShear_zoomed')
		Plot_stress_strain(-np.array(p),q,ce,pe,ee,'stress_strain_curve_noDamp_CompressionShear')


	# Jin routine model 

	if False:
		odb  = 'E-12.odb'
		path = 'G:/01_Forschung/01_Creep/00_UMAT/UMAT_CLAUDE/Comparing_Method/Comparison_Shengli/'
		p,q,t  = Extract_pq_ABQ_noUMAT(odb, path )
		ceeq = Extract_Strain_ABQ_CreepRoutine(odb, path)
		Plot_Pq(-np.array(p),q, 'dt='+str(dt))







delta_t = np.array([1,5,6,7,10,11,12,13,14])
dtt = 10.**-delta_t

if False:

	CPU_time = np.array([703,687,707,735,712,741,853,1620,8640])/60.

	cm = 1 / 2.54
	fig, axes = plt.subplots(figsize=(12*cm, 6*cm))

	axes.plot(dtt, CPU_time, '+-k')
	axes.set_ylabel(r"CPU time (min)")
	axes.set_yscale('log')
	axes.set_xscale('log')

	axes.set_xlabel(r"$\varepsilon_\mathrm{sub}$")
	axes.set_xlabel(r"$\varepsilon_0^\mathrm{cr}$")
	axes.set_yticks((10,100))
	axes.set_xlim(1e-15,2*10**-1)

	plt.tight_layout()

	plt.savefig(path+'cpu_time.png', dpi=300)
	
	plt.close("all")