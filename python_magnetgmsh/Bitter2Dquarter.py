import os
import sys

import yaml
import gmsh
from python_magnetgeo.Bitter import Bitter
from python_magnetgeo.Shape2D import Shape2D

from .utils.lists import flatten
from .mesh.axi import get_allowed_algo, MeshAlgo2D


# TieRod
def create_shape(x: float, y: float, shape: Shape2D):
    """
    create Shape
    """

    _cl = None
    print(f"shape: {shape}, name={shape.name}, x={x}, y={y}", flush=True)
    if shape.name.startswith("circle"):
        r = float(shape.name.split("-")[-1].replace("mm", "")) / 2.0
        curv = gmsh.model.occ.addCircle(x, y, 0, r)  # , angle1=, angle2 = )
        _cl = gmsh.model.occ.addCurveLoop([curv])
        print(f"create a circle: C({x},{y}), r={r}, _cl={_cl}, curv={curv}")
    else:
        curv = []
        pts = []
        for i, pt in enumerate(shape.pts):
            print(f"pt[{i}]: {pt}")
            pts.append(gmsh.model.occ.addPoint(x + pt[0], y + pt[1], 0))
            print(f"addPoint({x + pt[0]}, {y+pt[1]})", flush=True)
            if i >= 1:
                curv.append(gmsh.model.occ.addLine(pts[i - 1], pts[i]))
        curv.append(gmsh.model.occ.addLine(pts[-1], pts[0]))
        _cl = gmsh.model.occ.addCurveLoop(curv)
        print(f"create_shape: _cl={_cl}, {curv}")
        curv.clear()
        pts.clear()
    return _cl


def gmsh2D_ids(Bitter: Bitter, AirData: tuple, debug: bool = False) -> tuple:
    """
    create gmsh 2D geometry
    """
    print("Bitter/gmsh2D_ids")

    from math import pi, cos, sin

    theta = pi / 2
    if Bitter.tierod:
        tierod = Bitter.tierod

    Origin = gmsh.model.occ.addPoint(0, 0, 0)

    # Bitter sector
    curv = []
    pt0_r0 = gmsh.model.occ.addPoint(Bitter.r[0] * cos(0), Bitter.r[0] * sin(0), 0)
    pt1_r0 = gmsh.model.occ.addPoint(
        Bitter.r[0] * cos(theta), Bitter.r[0] * sin(theta), 0
    )
    rint_id = gmsh.model.occ.addCircleArc(pt0_r0, Origin, pt1_r0)
    print(f"rint_id={rint_id}")
    curv.append(rint_id)

    pt0_r1 = gmsh.model.occ.addPoint(Bitter.r[1] * cos(0), Bitter.r[1] * sin(0), 0)
    pt1_r1 = gmsh.model.occ.addPoint(
        Bitter.r[1] * cos(theta), Bitter.r[1] * sin(theta), 0
    )
    curv.append(gmsh.model.occ.addLine(pt0_r0, pt0_r1))
    rext_id = gmsh.model.occ.addCircleArc(pt0_r1, Origin, pt1_r1)
    print(f"rext_id={rext_id}")
    curv.append(rext_id)
    curv.append(gmsh.model.occ.addLine(pt1_r1, pt1_r0))
    cl = gmsh.model.occ.addCurveLoop(curv)
    sector = gmsh.model.occ.addPlaneSurface([cl])
    print(f"Bitter sector: {cl} = {curv}")
    del curv

    holes = []
    names = []
    tierods = []
    tnames = []
    if Bitter.tierod:
        ntierod = tierod.n

        theta_t = 2 * pi / float(ntierod)
        angle = pi / float(ntierod)
        tierod_thetas = [angle]

        # create Shape for slit
        _ltierod = create_shape(tierod.r, 0, tierod.shape)
        print(f"_ltierod: {_ltierod}", flush=True)
        tierod_id = gmsh.model.occ.addPlaneSurface([_ltierod])
        gmsh.model.occ.rotate([(2, tierod_id)], 0, 0, 0, 0, 0, 1, angle)
        print(f"tierod[0]: rotate {angle} init")
        tierods.append(tierod_id)
        tnames.append("tierod_0")

        for n in range(1, ntierod):
            if (
                n * theta_t + angle <= theta
                or n * theta_t + angle >= 2 * pi  # - theta / 2.0
            ):
                res = gmsh.model.occ.copy([(2, tierod_id)])
                # print(f"res={res}")
                _id = res[0][1]
                gmsh.model.occ.rotate([(2, _id)], 0, 0, 0, 0, 0, 1, n * theta_t)
                tierods.append(_id)
                tnames.append(f"tierod_{n}")
                tierod_thetas.append(n * theta_t + angle)

        print(tierod_thetas)

        # gmsh.model.occ.rotate([(2, tierod_id)], 0, 0, 0, 0, 0, 1, theta/2.0)
        print(f"tierod_id: {tierod_id}", flush=True)

    # CoolingSlits
    if Bitter.coolingslits:
        for j, slit in enumerate(Bitter.coolingslits):
            _names = []
            print(
                f"slit[{j+1}]: nslits={slit.n}, r={slit.r}",  # tierod={tierod.r* cos(theta/2.0)}",
                flush=True,
            )
            nslits = slit.n

            if Bitter.tierod and slit.r == tierod.r:
                nslits += tierod.n

            theta_s = 2 * pi / float(nslits)
            angle = slit.angle * pi / 180.0

            # create Shape for slit
            _lc = create_shape(x=slit.r, y=0, shape=slit.shape)
            slit_id = gmsh.model.occ.addPlaneSurface([_lc])
            if angle != 0:
                gmsh.model.occ.rotate([(2, slit_id)], 0, 0, 0, 0, 0, 1, angle)
                print(f"slit[{j+1}][0]: rotate {angle} init")

            if Bitter.tierod and slit.r == tierod.r and angle == theta / 2.0:
                print(f"angle={angle}, tierod_thetas={tierod_thetas}")
                print("skip slit")
            else:
                holes.append(slit_id)
                _names.append(f"slit{j+1}_0")

            for n in range(1, nslits):
                if (
                    n * theta_s + angle <= theta
                    or n * theta_s + angle >= 2 * pi  # - theta / 2.0
                ):
                    if (
                        Bitter.tierod
                        and slit.r == tierod.r
                        and n * theta_s + angle in tierod_thetas
                    ):
                        print("skip slit")
                    else:
                        res = gmsh.model.occ.copy([(2, slit_id)])
                        # print(f"res={res}")
                        _id = res[0][1]
                        gmsh.model.occ.rotate([(2, _id)], 0, 0, 0, 0, 0, 1, n * theta_s)
                        holes.append(_id)
                        _names.append(f"slit{j+1}_{n}")

            names.append(_names)

            # if slit.r == tierod.r and angle == theta / 2.0:
            #     print(f"remove slit{j+1}_0: {slit_id}")
            #     gmsh.model.occ.remove([(2, slit_id)], recursive=True)
            #     gmsh.model.occ.synchronize()

    if Bitter.tierod or Bitter.coolingslits:
        slit_tierod = flatten(holes) + flatten(tierods)
        name_slit_tierod = flatten(names) + flatten(tnames)
        print(f"slits+tierods={slit_tierod}")
        print(f"names={name_slit_tierod}")
        cad = gmsh.model.occ.cut(
            [(2, sector)], [(2, _id) for _id in slit_tierod], removeTool=False
        )
    else:
        cad = sector
    gmsh.model.occ.synchronize()

    if holes or tierods:
        print(f"cad: {cad}")
        print(f"cad[1]: {cad[1]}")
        print(f"cad[1][0]: {cad[1][0]}")
        print(f"PhysicalGroup[{Bitter.name}]: ids={[cad[0][0][1]]}", flush=True)
        ps = gmsh.model.addPhysicalGroup(2, [cad[0][0][1]])
        gmsh.model.setPhysicalName(2, ps, Bitter.name)
        print(f"PhysicalGroup[{Bitter.name}]: ids={[cad[0][0][1]]}", flush=True)
    else:
        print(f"cad: {cad}")
        ps = gmsh.model.addPhysicalGroup(2, [cad])
        gmsh.model.setPhysicalName(2, ps, Bitter.name)
        print(f"PhysicalGroup[{Bitter.name}]: ids={[cad]} (sector only)", flush=True)

    # get BCs ids
    # use
    # gmsh/model/occ/getBoundingBox
    # gmsh/model/occ/getEntitiesInBoundingBox
    def create_bcgroup(shape: int, subshape: int, name: str):
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.occ.getBoundingBox(2, subshape)
        # print(f"boundingbox[{name}]: {[xmin, ymin, zmin, xmax, ymax, zmax]}")
        if abs(zmin - zmax) >= 1.0e-6:
            raise RuntimeError(
                f"create_bcgroup({name}): subshape is expected to be in OxOy plane (zmin={zmin}, Zmax={zmax})"
            )

        interface = gmsh.model.occ.getEntitiesInBoundingBox(
            xmin, ymin, zmin, xmax, ymax, zmax, dim=1
        )

        ids = [tag[1] for tag in interface]
        print(f"interface[{name}]: ids={ids}, interface={interface}")
        ps = gmsh.model.addPhysicalGroup(1, ids)
        gmsh.model.setPhysicalName(1, ps, name)

        return ids

    # gmsh.model.occ.removeAllDuplicates()
    # print(tierod_id)
    # tierod_ids = create_bcgroup(cad[0][0][1], tierod_id, "tierod")

    slit_names = flatten(names)
    print(f"slit_names: {len(slit_names)} names, {len(holes)} slits")
    hole_ids = []
    for i in range(0, len(holes)):
        _ids = create_bcgroup(cad[0][0][1], holes[i], slit_names[i])
        hole_ids.append(_ids)

    print(f"tierod_names: {len(tnames)} names, {len(tierods)} tierods")
    tierod_ids = []
    for i in range(0, len(tierods)):
        _ids = create_bcgroup(cad[0][0][1], tierods[i], tnames[i])
        tierod_ids.append(_ids)

    for hole in holes:
        gmsh.model.occ.remove([(2, hole)], recursive=True)

    for tierod in tierods:
        gmsh.model.occ.remove([(2, tierod)], recursive=True)

    gmsh.model.occ.synchronize()

    # for tierod in tierod_ids:
    #     print(tierod)
    #     gmsh.model.occ.remove([(2, tierod)], recursive=True)
    # gmsh.model.occ.synchronize()

    # TODO physical for Rint, Rext, V0, V1
    eps = 0.1
    xmin = Bitter.r[0] * cos(theta) - eps
    ymin = -eps
    xmax = Bitter.r[0] + eps
    ymax = Bitter.r[0] * sin(theta) + eps
    candidate_rint = gmsh.model.getEntitiesInBoundingBox(
        xmin,
        ymin,
        -eps,
        xmax,
        ymax,
        eps,
        1,
    )
    candidate_rint_ids = [tag[1] for tag in candidate_rint]
    print(
        f"candidate_rint={candidate_rint} xmin={xmin} ymin={ymin} xmax={xmax}, ymax={ymax}"
    )
    xmin = Bitter.r[1] * cos(theta) - eps
    ymin = -eps
    xmax = Bitter.r[1] + eps
    ymax = Bitter.r[1] * sin(theta) + eps
    candidate_rext = gmsh.model.getEntitiesInBoundingBox(
        xmin,
        ymin,
        -eps,
        xmax,
        ymax,
        eps,
        1,
    )
    candidate_rext_ids = [tag[1] for tag in candidate_rext]
    print(
        f"candidate_rext={candidate_rext} xmin={xmin} ymin={ymin} xmax={xmax}, ymax={ymax}"
    )
    xmin = Bitter.r[0] * cos(0) - eps
    ymin = Bitter.r[1] * sin(0) - eps
    xmax = Bitter.r[1] * cos(0) + eps
    ymax = Bitter.r[0] * sin(0) + eps
    candidate_V0 = gmsh.model.getEntitiesInBoundingBox(
        xmin,
        ymin,
        -eps,
        xmax,
        ymax,
        eps,
        1,
    )
    candidate_V0_ids = [tag[1] for tag in candidate_V0]
    print(
        f"candidate_V0={candidate_V0} xmin={xmin} ymin={ymin} xmax={xmax}, ymax={ymax}"
    )
    xmin = Bitter.r[0] * cos(theta) - eps
    ymin = Bitter.r[0] * sin(theta) - eps
    xmax = Bitter.r[1] * cos(theta) + eps
    ymax = Bitter.r[1] * sin(theta) + eps
    candidate_V1 = gmsh.model.getEntitiesInBoundingBox(
        xmin,
        ymin,
        -eps,
        xmax,
        ymax,
        eps,
        1,
    )
    candidate_V1_ids = [tag[1] for tag in candidate_V1]
    print(
        f"candidate_V1={candidate_V1} xmin={xmin} ymin={ymin} xmax={xmax}, ymax={ymax}"
    )

    _ids = flatten(hole_ids) + flatten(tierod_ids)
    print(_ids)
    vEntities = gmsh.model.getEntities(1)
    rint_ids = []
    rext_ids = []
    V0_ids = []
    V1_ids = []
    for i, entity in enumerate(vEntities):
        tag = entity[1]
        if entity[1] not in _ids:
            gtype = gmsh.model.getType(entity[0], entity[1])
            print(f"Line[{i}]: id={tag}, type={gtype}")
            if gtype == "Circle":
                if tag in candidate_rint_ids:
                    rint_ids.append(tag)
                if tag in candidate_rext_ids and tag not in candidate_rint_ids:
                    rext_ids.append(tag)
            elif gtype == "Line":
                if tag in candidate_V0_ids:
                    V0_ids.append(tag)
                if tag in candidate_V1_ids:
                    V1_ids.append(tag)

    ps = gmsh.model.addPhysicalGroup(1, rint_ids)
    gmsh.model.setPhysicalName(1, ps, "slit0")
    ps = gmsh.model.addPhysicalGroup(1, rext_ids)
    if Bitter.coolingslits:
        gmsh.model.setPhysicalName(1, ps, f"slit{len(Bitter.coolingslits)+1}")
    else:
        gmsh.model.setPhysicalName(1, ps, "slit1")
    ps = gmsh.model.addPhysicalGroup(1, V0_ids)
    gmsh.model.setPhysicalName(1, ps, "V0")
    ps = gmsh.model.addPhysicalGroup(1, V1_ids)
    gmsh.model.setPhysicalName(1, ps, "V1")

    if isinstance(cad, int):
        return (
            [cad],
            (rint_ids + rext_ids + V0_ids + V1_ids),
            flatten(hole_ids),
            flatten(tierod_ids),
        )
    else:
        return (
            [cad[1][0][0][1]],
            (rint_ids + rext_ids + V0_ids + V1_ids),
            flatten(hole_ids),
            flatten(tierod_ids),
        )


def main():
    import argparse

    """Console script for python_magnetgeo."""
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "filename", help="name of the model to be loaded (yaml file)", type=str
    )
    parser.add_argument("--wd", help="set a working directory", type=str, default="")
    parser.add_argument("--mesh", help="activate mesh", action="store_true")
    parser.add_argument(
        "--algo2d",
        help="select an algorithm for 2d mesh",
        type=str,
        choices=get_allowed_algo(),
        default="Delaunay",
    )
    parser.add_argument(
        "--scaling", help="scale to m (default unit is mm)", action="store_true"
    )
    parser.add_argument("--lc", help="specify mesh size", type=float, default="10")
    parser.add_argument("--show", help="display gmsh windows", action="store_true")
    parser.add_argument("--verbose", help="activate debug mode", action="store_true")
    parser.add_argument("--debug", help="activate debug mode", action="store_true")
    args = parser.parse_args()
    print(f"Arguments: {args}")

    cwd = os.getcwd()
    if args.wd:
        os.chdir(args.wd)

    with open(args.filename, "r") as f:
        Object = yaml.load(f, Loader=yaml.FullLoader)
        print(f"Object={Object}, type={type(Object)}")

    ncoolingslits = 0
    if Object.coolingslits:
        ncoolingslits += len(Object.coolingslits)
    print(f"ncoolingslits={ncoolingslits}")

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("General.Verbosity", 0)
    if args.debug or args.verbose:
        gmsh.option.setNumber("General.Verbosity", 5)
    gmsh.model.add(args.filename)
    gmsh.logger.start()

    AirData = ()
    result = gmsh2D_ids(Object, AirData, args.debug)
    print(f"gmsh2D_ids={result}")
    # Todo: create BCs

    if args.mesh:
        print(f"create Axi Gmsh mesh ({args.algo2d})")
        gmsh.option.setNumber("Mesh.Algorithm", MeshAlgo2D[args.algo2d])

        # scaling
        unit = 1
        if args.scaling:
            unit = 0.001
            gmsh.option.setNumber("Geometry.OCCScaling", unit)

        gmsh.model.mesh.setSize(gmsh.model.getEntities(0), args.lc)

        # adapt mesh
        dfields = []
        nfield = 0

        # Rint/Rext
        for bc in result[1]:
            gmsh.model.mesh.field.add("Distance", nfield)
            gmsh.model.mesh.field.setNumbers(nfield, "CurvesList", [bc])
            gmsh.model.mesh.field.setNumbers(nfield, "Sampling", [100])
            nfield += 1

            # Field 2: Threshold that dictates the mesh size of the background field
            print(f"Field[Thresold] for cooling holes: {nfield}")
            gmsh.model.mesh.field.add("Threshold", nfield)
            gmsh.model.mesh.field.setNumber(nfield, "IField", nfield - 1)
            gmsh.model.mesh.field.setNumber(nfield, "LcMin", args.lc / 10.0 * unit)
            gmsh.model.mesh.field.setNumber(nfield, "LcMax", args.lc * unit)
            gmsh.model.mesh.field.setNumber(nfield, "DistMin", 10 * unit)
            gmsh.model.mesh.field.setNumber(nfield, "DistMax", 12.5 * unit)
            gmsh.model.mesh.field.setNumber(nfield, "StopAtDistMax", True)
            dfields.append(nfield)
            nfield += 1

        # Tierod
        if Object.tierod:
            r_tierod = min(0.5, Object.tierod.r)
            for bc in result[3]:
                gmsh.model.mesh.field.add("Distance", nfield)
                gmsh.model.mesh.field.setNumbers(nfield, "CurvesList", [bc])
                gmsh.model.mesh.field.setNumbers(nfield, "Sampling", [100])
                nfield += 1

                # Field 2: Threshold that dictates the mesh size of the background field
                print(f"Field[Thresold] for cooling holes: {nfield}")
                gmsh.model.mesh.field.add("Threshold", nfield)
                gmsh.model.mesh.field.setNumber(nfield, "IField", nfield - 1)
                gmsh.model.mesh.field.setNumber(nfield, "LcMin", args.lc / 20.0 * unit)
                gmsh.model.mesh.field.setNumber(nfield, "LcMax", args.lc * unit)
                gmsh.model.mesh.field.setNumber(nfield, "DistMin", r_tierod * unit)
                gmsh.model.mesh.field.setNumber(
                    nfield, "DistMax", 1.1 * r_tierod * unit
                )
                gmsh.model.mesh.field.setNumber(nfield, "StopAtDistMax", True)
                dfields.append(nfield)
                nfield += 1

        # Cooling holes
        lst = []
        r_holes = 2
        if Object.coolingslits:
            lst = [slit.dh for slit in Object.coolingslits]
            r_holes = min(2, sum(lst) / len(lst))
        for bc in result[2]:
            gmsh.model.mesh.field.add("Distance", nfield)
            gmsh.model.mesh.field.setNumbers(nfield, "CurvesList", [bc])
            gmsh.model.mesh.field.setNumbers(nfield, "Sampling", [100])
            nfield += 1

            # Field 2: Threshold that dictates the mesh size of the background field
            print(f"Field[Thresold] for cooling holes: {nfield}")
            gmsh.model.mesh.field.add("Threshold", nfield)
            gmsh.model.mesh.field.setNumber(nfield, "IField", nfield - 1)
            gmsh.model.mesh.field.setNumber(nfield, "LcMin", args.lc / 40.0 * unit)
            gmsh.model.mesh.field.setNumber(nfield, "LcMax", args.lc * unit)
            gmsh.model.mesh.field.setNumber(nfield, "DistMin", r_holes * unit)
            gmsh.model.mesh.field.setNumber(nfield, "DistMax", 1.02 * r_holes * unit)
            gmsh.model.mesh.field.setNumber(nfield, "StopAtDistMax", True)
            dfields.append(nfield)
            nfield += 1

        # Let's use the minimum of all the fields as the mesh size field:
        gmsh.model.mesh.field.add("Min", nfield)
        gmsh.model.mesh.field.setNumbers(nfield, "FieldsList", dfields)  # dfields
        print(f"Field[Min] = {nfield}, Min=[{dfields[0]},...,{dfields[-1]}]")

        print(f"Apply background mesh {nfield}")
        gmsh.model.mesh.field.setAsBackgroundMesh(nfield)

        gmsh.model.mesh.generate(2)

        # TODO create 2D mesh
        gmsh.option.setNumber("Mesh.SaveAll", 1)
        meshfilename = args.filename.replace(".yaml", "-2D")
        gmsh.write(meshfilename + ".msh")

    gmsh.model.occ.synchronize()
    if args.show:
        gmsh.fltk.run()
    gmsh.finalize()


if __name__ == "__main__":
    sys.exit(main())
