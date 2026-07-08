# -*- coding: utf-8 -*-

from odbAccess import openOdb
from abaqusConstants import INTEGRATION_POINT
import os
import csv
import math



#abaqus cae noGUI=C:\Users\p2515497\Documents\ghk-creep\ghk-creep\inp_schachner.py

# ============================================================
# USER INPUTS
# ============================================================

ODB_PATH = r"G:\01_Forschung\01_Creep\00_UMAT\UMAT_CLAUDE\Schachter_3D\G1_schachter.odb"

STEP_NAME = "Load1"     # mets None pour prendre le dernier step
FRAME_ID = -1            # -1 = dernière frame
OUT_DAT = r"C:\Users\p2515497\Documents\ghk-creep\ghk-creep\sdv3_sdv4_ip_ecr0_8.dat"
OUT_DAT_Path = r"C:\Users\p2515497\Documents\ghk-creep\ghk-creep\sdv3_sdv4_path_ecr0_8.dat"

# ============================================================
# FUNCTIONS
# ============================================================

def field_by_ip_key(field):
    """
    Build dictionary:
        key = (instanceName, elementLabel, integrationPoint)
        value = scalar data
    """
    d = {}

    for v in field.values:
        inst_name = v.instance.name
        elem_label = v.elementLabel
        ip = v.integrationPoint
        val = float(v.data)

        d[(inst_name, elem_label, ip)] = val

    return d


def coord_by_ip_from_coord_field(frame):
    """
    Try to read integration point coordinates directly from COORD field.
    This is the best method, but COORD must be present in the ODB.
    """
    if "COORD" not in frame.fieldOutputs.keys():
        return None

    coord_field = frame.fieldOutputs["COORD"].getSubset(position=INTEGRATION_POINT)

    d = {}
    for v in coord_field.values:
        inst_name = v.instance.name
        elem_label = v.elementLabel
        ip = v.integrationPoint

        data = tuple(float(a) for a in v.data)

        if len(data) == 2:
            data = (data[0], data[1], 0.0)

        d[(inst_name, elem_label, ip)] = data

    return d


def build_node_coord_cache(instance):
    d = {}
    for n in instance.nodes:
        c = tuple(float(a) for a in n.coordinates)
        if len(c) == 2:
            c = (c[0], c[1], 0.0)
        d[n.label] = c
    return d


def interp_coord(node_coords, N):
    x = 0.0
    y = 0.0
    z = 0.0

    for a, Na in enumerate(N):
        x += Na * node_coords[a][0]
        y += Na * node_coords[a][1]
        z += Na * node_coords[a][2]

    return (x, y, z)


def shape_C3D10(L1, L2, L3, L4):
    """
    Quadratic tetrahedral C3D10 shape functions.

    Abaqus C3D10 node ordering:
      1,2,3,4 = corner nodes
      5 = edge 1-2
      6 = edge 2-3
      7 = edge 3-1
      8 = edge 1-4
      9 = edge 2-4
      10 = edge 3-4
    """
    return [
        L1 * (2.0 * L1 - 1.0),
        L2 * (2.0 * L2 - 1.0),
        L3 * (2.0 * L3 - 1.0),
        L4 * (2.0 * L4 - 1.0),
        4.0 * L1 * L2,
        4.0 * L2 * L3,
        4.0 * L3 * L1,
        4.0 * L1 * L4,
        4.0 * L2 * L4,
        4.0 * L3 * L4,
    ]


def coord_by_ip_reconstructed(odb, sdv_keys):
    """
    Reconstruct integration point coordinates for C3D10 elements.

    This is used only if COORD is not available in the ODB.
    """
    asm = odb.rootAssembly

    node_cache = {}
    elem_cache = {}

    # 4 integration points of tetrahedron
    a_tet = 0.5854101966249685
    b_tet = 0.1381966011250105

    gps_C3D10 = [
        (a_tet, b_tet, b_tet, b_tet),
        (b_tet, a_tet, b_tet, b_tet),
        (b_tet, b_tet, a_tet, b_tet),
        (b_tet, b_tet, b_tet, a_tet),
    ]

    out = {}

    for key in sdv_keys:
        inst_name, elem_label, ip = key

        if inst_name not in asm.instances.keys():
            out[key] = ("", "", "")
            continue

        inst = asm.instances[inst_name]

        if inst_name not in node_cache:
            node_cache[inst_name] = build_node_coord_cache(inst)

        if (inst_name, elem_label) not in elem_cache:
            try:
                elem_cache[(inst_name, elem_label)] = inst.getElementFromLabel(elem_label)
            except:
                out[key] = ("", "", "")
                continue

        elem = elem_cache[(inst_name, elem_label)]
        etype = elem.type.upper()
        conn = elem.connectivity

        nodes = [node_cache[inst_name][lab] for lab in conn]

        if etype in ["C3D10", "C3D10M"]:
            if ip >= 1 and ip <= 4:
                L1, L2, L3, L4 = gps_C3D10[ip - 1]
                N = shape_C3D10(L1, L2, L3, L4)
                out[key] = interp_coord(nodes, N)
            else:
                out[key] = ("", "", "")
        else:
            out[key] = ("", "", "")

    return out

from odbAccess import openOdb
import math


# ============================================================
# Small vector tools, Abaqus-safe
# ============================================================

def vec_sub(a, b):
    return tuple(a[i] - b[i] for i in range(len(a)))


def vec_add(a, b):
    return tuple(a[i] + b[i] for i in range(len(a)))


def vec_mul(alpha, a):
    return tuple(alpha * a[i] for i in range(len(a)))


def vec_dot(a, b):
    s = 0.0
    for i in range(len(a)):
        s += a[i] * b[i]
    return s


def vec_norm(a):
    s = 0.0
    for i in range(len(a)):
        s += a[i] * a[i]
    return math.sqrt(s)


def dist(a, b):
    return vec_norm(vec_sub(a, b))


# ============================================================
# SDV tools
# ============================================================

def get_sdv_field(frame, sdv_number):
    """
    Works for:
        SDV3, SDV4
    or for:
        SDV with components.
    """

    key = "SDV{}".format(sdv_number)

    if key in frame.fieldOutputs.keys():
        return frame.fieldOutputs[key], None

    if "SDV" in frame.fieldOutputs.keys():
        return frame.fieldOutputs["SDV"], sdv_number - 1

    raise RuntimeError("Could not find SDV{} or SDV field.".format(sdv_number))


def get_sdv_value(value, component_index):
    if component_index is None:
        return float(value.data)
    else:
        return float(value.data[component_index])


# ============================================================
# Shape functions
# ============================================================

def centroid(coords):
    ndim = len(coords[0])
    c = [0.0] * ndim

    for x in coords:
        for i in range(ndim):
            c[i] += x[i]

    for i in range(ndim):
        c[i] /= float(len(coords))

    return tuple(c)


def q4_shape_functions(xi, eta):
    return [
        0.25 * (1.0 - xi) * (1.0 - eta),
        0.25 * (1.0 + xi) * (1.0 - eta),
        0.25 * (1.0 + xi) * (1.0 + eta),
        0.25 * (1.0 - xi) * (1.0 + eta),
    ]


def c3d8_shape_functions(xi, eta, zeta):
    return [
        0.125 * (1.0 - xi) * (1.0 - eta) * (1.0 - zeta),
        0.125 * (1.0 + xi) * (1.0 - eta) * (1.0 - zeta),
        0.125 * (1.0 + xi) * (1.0 + eta) * (1.0 - zeta),
        0.125 * (1.0 - xi) * (1.0 + eta) * (1.0 - zeta),
        0.125 * (1.0 - xi) * (1.0 - eta) * (1.0 + zeta),
        0.125 * (1.0 + xi) * (1.0 - eta) * (1.0 + zeta),
        0.125 * (1.0 + xi) * (1.0 + eta) * (1.0 + zeta),
        0.125 * (1.0 - xi) * (1.0 + eta) * (1.0 + zeta),
    ]


def interpolate_from_shape_functions(N, node_coords):
    ndim = len(node_coords[0])
    x = [0.0] * ndim

    for a in range(len(N)):
        for i in range(ndim):
            x[i] += N[a] * node_coords[a][i]

    return tuple(x)


def compute_element_ip_coords(node_coords, n_ip):
    """
    Reconstruct integration point coordinates.

    Supported:
        Q4 full integration: 4 IP
        C3D8 full integration: 8 IP
        reduced integration: 1 IP -> centroid
    """

    if n_ip == 1:
        return [centroid(node_coords)]

    g = 1.0 / math.sqrt(3.0)

    # 2D Q4, 4 integration points
    if len(node_coords) == 4 and n_ip == 4:
        gps = [
            (-g, -g),
            ( g, -g),
            ( g,  g),
            (-g,  g),
        ]

        coords_ip = []

        for xi, eta in gps:
            N = q4_shape_functions(xi, eta)
            coords_ip.append(interpolate_from_shape_functions(N, node_coords))

        return coords_ip

    # 3D C3D8, 8 integration points
    if len(node_coords) == 8 and n_ip == 8:
        gps = [
            (-g, -g, -g),
            ( g, -g, -g),
            ( g,  g, -g),
            (-g,  g, -g),
            (-g, -g,  g),
            ( g, -g,  g),
            ( g,  g,  g),
            (-g,  g,  g),
        ]

        coords_ip = []

        for xi, eta, zeta in gps:
            N = c3d8_shape_functions(xi, eta, zeta)
            coords_ip.append(interpolate_from_shape_functions(N, node_coords))

        return coords_ip

    # Fallback: use centroid for every IP
    c = centroid(node_coords)
    return [c for i in range(n_ip)]

from odbAccess import openOdb
from abaqus import session
from abaqusConstants import *


def extract_sdv_from_abaqus_path(
        odb_path,
        point_1,
        point_2,
        step_name,
        sdva,sdvb,
        frame_id=-1,
        n_points=100,
        output_dat="path_sdv3_sdv4.dat",
        ):
    """
    Extract SDV3 and SDV4 along a straight Abaqus path,
    using Abaqus' own Data from Path algorithm.

    This should match:
        Tools -> Path -> Create
        Tools -> XY Data -> Create -> Path
        Shape: Undeformed
        Number of points: 100
        Average value threshold: 50%
    """

    # ------------------------------------------------------------
    # Open ODB in CAE session
    # ------------------------------------------------------------
    odb = session.openOdb(name=odb_path)

    vp_name = "Viewport: 1"

    if vp_name not in session.viewports.keys():
        session.Viewport(name=vp_name)

    vp = session.viewports[vp_name]
    vp.setValues(displayedObject=odb)

    # ------------------------------------------------------------
    # Set step and frame
    # ------------------------------------------------------------
    step_names = odb.steps.keys()

    if step_name not in step_names:
        odb.close()
        raise RuntimeError("Step '{}' not found. Available steps: {}".format(
            step_name, step_names
        ))

    step_index = step_names.index(step_name)

    vp.odbDisplay.setFrame(step=step_index, frame=frame_id)

    # ------------------------------------------------------------
    # Create path from two points
    # ------------------------------------------------------------
    p1 = tuple(float(x) for x in point_1)
    p2 = tuple(float(x) for x in point_2)

    if len(p1) == 2:
        p1 = (p1[0], p1[1], 0.0)

    if len(p2) == 2:
        p2 = (p2[0], p2[1], 0.0)

    path_name = "MY_SCRIPT_PATH"

    if path_name in session.paths.keys():
        del session.paths[path_name]

    path_obj = session.Path(
        name=path_name,
        type=POINT_LIST,
        expression=(p1, p2)
    )

    # Abaqus numIntervals means number of intervals,
    # so 100 points = 99 intervals.
    num_intervals = n_points - 1

    # ------------------------------------------------------------
    # Extract SDV3
    # ------------------------------------------------------------
    vp.odbDisplay.setPrimaryVariable(
        variableLabel=sdva,
        outputPosition=INTEGRATION_POINT
    )

    xy_sdv3 = session.XYDataFromPath(
        name=sdva+"_path",
        path=path_obj,
        includeIntersections=False,
        pathStyle=UNIFORM_SPACING,
        numIntervals=num_intervals,
        shape=UNDEFORMED,
        labelType=TRUE_DISTANCE,
        removeDuplicateXYPairs=True,
        includeAllElements=False
    )

    # ------------------------------------------------------------
    # Extract SDV4
    # ------------------------------------------------------------
    vp.odbDisplay.setPrimaryVariable(
        variableLabel=sdvb,
        outputPosition=INTEGRATION_POINT
    )

    xy_sdv4 = session.XYDataFromPath(
        name=sdvb+"_path",
        path=path_obj,
        includeIntersections=False,
        pathStyle=UNIFORM_SPACING,
        numIntervals=num_intervals,
        shape=UNDEFORMED,
        labelType=TRUE_DISTANCE,
        removeDuplicateXYPairs=True,
        includeAllElements=False
    )

    # ------------------------------------------------------------
    # Write .dat file
    # ------------------------------------------------------------
    data3 = list(xy_sdv3.data)
    data4 = list(xy_sdv4.data)

    n = min(len(data3), len(data4))

    f = open(output_dat, "w")
    f.write("# distance SDV3 SDV4\n")

    for i in range(n):
        s3, v3 = data3[i]
        s4, v4 = data4[i]

        # They should be the same distance, but to be safe:
        s = 0.5 * (s3 + s4)

        f.write("{:.8e} {:.8e} {:.8e}\n".format(s, v3, v4))

    f.close()

    odb.close()

    return data3, data4


def Extract_force_time(odb_path, step_name, set_extraction, direction_force, out_file):

    odb         = openOdb(odb_path, readOnly=True)

    step        = odb.steps[step_name]

    if set_extraction in odb.rootAssembly.nodeSets.keys():
        node_set = odb.rootAssembly.nodeSets[set_extraction]
    else:
        odb.close()
        raise RuntimeError("Node set '{} not found in rootAssembly".format(set_extraction))

    f=open(out_file,"w")
    f.write(" time \t RF \n")

    for frame in step.frames:
        time    = frame.frameValue

        if "RF" not in frame.fieldOutputs.keys():
            odb.close()
            raise RuntimeError("RF field output not found. Please request in the .inp")

        rf_field    = frame.fieldOutputs["RF"]
        rf_subject  = rf_field.getSubset(region=node_set)

        total_force = 0.

        for value in rf_subject.values:
            total_force += value.data[direction_force]
        
        f.write("{:.4e} \t {:.2e} \n".format(time, total_force))
    
    f.close()
    odb.close()





# ============================================================
# MAIN
# ============================================================

def main():

    print("Opening ODB:")
    print(ODB_PATH)

    odb = openOdb(ODB_PATH, readOnly=True)

    try:
        if STEP_NAME is None:
            step_name = list(odb.steps.keys())[-1]
        else:
            step_name = STEP_NAME

        if step_name not in odb.steps.keys():
            raise RuntimeError("Step '%s' not found. Available steps: %s"
                               % (step_name, list(odb.steps.keys())))

        step = odb.steps[step_name]
        frame = step.frames[FRAME_ID]

        print("Step:  ", step.name)
        print("Frame: ", FRAME_ID)
        print("Time:  ", frame.frameValue)

        if "SDV3" not in frame.fieldOutputs.keys():
            raise RuntimeError("SDV3 not found in field outputs.")

        if "SDV4" not in frame.fieldOutputs.keys():
            raise RuntimeError("SDV4 not found in field outputs.")

        sdv3_field = frame.fieldOutputs["SDV3"].getSubset(position=INTEGRATION_POINT)
        sdv4_field = frame.fieldOutputs["SDV4"].getSubset(position=INTEGRATION_POINT)

        sdv3 = field_by_ip_key(sdv3_field)
        sdv4 = field_by_ip_key(sdv4_field)

        keys = sorted(set(sdv3.keys()) & set(sdv4.keys()))

        print("Number of integration point values:", len(keys))

        coords = coord_by_ip_from_coord_field(frame)

        if coords is None:
            print("COORD at integration points not found.")
            print("Reconstructing IP coordinates for C3D10/C3D10M.")
            coords = coord_by_ip_reconstructed(odb, keys)
        else:
            print("Using COORD field at integration points.")

        with open(OUT_DAT, "w") as f:

            f.write("# Integration point output from Abaqus ODB\n")
            f.write("# ODB:   %s\n" % ODB_PATH)
            f.write("# Step:  %s\n" % step.name)
            f.write("# Frame: %s\n" % str(FRAME_ID))
            f.write("# Time:  %.12e\n" % frame.frameValue)
            f.write("#\n")
            f.write("# Columns:\n")
            f.write("# instance  elementLabel  integrationPoint  x  y  z  SDV3  SDV4\n")
            f.write("#\n")

            for key in keys:
                inst_name, elem_label, ip = key
                x, y, z = coords.get(key, ("", "", ""))

                if x == "":
                    f.write("%s %d %d NaN NaN NaN %.12e %.12e\n" %
                            (inst_name, elem_label, ip, sdv3[key], sdv4[key]))
                else:
                    f.write("%s %d %d %.12e %.12e %.12e %.12e %.12e\n" %
                            (inst_name, elem_label, ip,
                            x, y, z,
                            sdv3[key], sdv4[key]))

        print("Done.")
        print("Written DAT:")
        print(OUT_DAT)

    finally:
        odb.close()

from abaqus import *
from abaqusConstants import *
from caeModules import *
from visualization import *
from driverUtils import executeOnCaeStartup

executeOnCaeStartup()

if __name__ == "__main__":
    if False:
        main()



    if False:
        point_1 = (-11.25,-15.0,           97.9807586669922)
        point_2 = (-11.25,15.0,            52.0192375183105)
        point_3 = (0.0,-15.0,97.9807586669922)
        point_4 = (0.0,15.0,52.0192375183105)


        OUT_DAT_Path1 = r"C:\Users\p2515497\Documents\ghk-creep\ghk-creep\sdv3_sdv4_path_surface_ecr0_8.dat"
        OUT_DAT_Path2 = r"C:\Users\p2515497\Documents\ghk-creep\ghk-creep\sdv3_sdv4_path_middle_ecr0_8.dat"


        if True:
            data_path = extract_sdv_from_abaqus_path(
                odb_path=ODB_PATH,
                point_1=point_1,
                point_2=point_2,
                step_name=STEP_NAME,
                sdva="SDV3",sdvb="SDV4",
                frame_id=-1,
                n_points=100,
                output_dat=OUT_DAT_Path1
            )
        

            data_path = extract_sdv_from_abaqus_path(
                odb_path=ODB_PATH,
                point_1=point_3,
                point_2=point_4,
                step_name=STEP_NAME,
                sdva="SDV3",sdvb="SDV4",
                frame_id=-1,
                n_points=100,
                output_dat=OUT_DAT_Path2
            )


        OUT_DAT_Path1 = r"C:\Users\p2515497\Documents\ghk-creep\ghk-creep\sdv1_sdv2_path_surface_ecr0_8.dat"
        OUT_DAT_Path2 = r"C:\Users\p2515497\Documents\ghk-creep\ghk-creep\sdv1_sdv2_path_middle_seed10.dat"
        if True:
            data_path = extract_sdv_from_abaqus_path(
                odb_path=ODB_PATH,
                point_1=point_1,
                point_2=point_2,
                step_name=STEP_NAME,
                sdva="SDV1",sdvb="SDV2",
                frame_id=-1,
                n_points=100,
                output_dat=OUT_DAT_Path1
            )
        

        data_path = extract_sdv_from_abaqus_path(
            odb_path=ODB_PATH,
            point_1=point_3,
            point_2=point_4,
            step_name=STEP_NAME,
            sdva="SDV1",sdvb="SDV2",
            frame_id=-1,
            n_points=100,
            output_dat=OUT_DAT_Path2
        )

    if True:
        out_fil = r"C:\Users\p2515497\Documents\ghk-creep\ghk-creep\Reac_force_seed1.dat"
        Extract_force_time(ODB_PATH,STEP_NAME, "TOP",2,out_fil)
