#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
Provides definition for Helix:

* Geom data: r, z
* Model Axi: definition of helical cut (provided from MagnetTools)
* Model 3D: actual 3D CAD
* Shape: definition of Shape eventually added to the helical cut
"""
from python_magnetgeo.Helix import Helix

import gmsh
from .mesh.bcs import create_bcs


def gmsh_ids(Helix: Helix, AirData: tuple, debug: bool = False) -> tuple:
    """
    create gmsh geometry
    """

    # TODO get axi model
    gmsh_ids = []
    x = Helix.r[0]
    dr = Helix.r[1] - Helix.r[0]
    y = -Helix.modelaxi.h

    if abs(y - Helix.z[0]) >= 0:
        _id = gmsh.model.occ.addRectangle(x, Helix.z[0], 0, dr, abs(y - Helix.z[0]))
        gmsh_ids.append(_id)

    for i, (n, pitch) in enumerate(zip(Helix.modelaxi.turns, Helix.modelaxi.pitch)):
        dz = n * pitch
        _id = gmsh.model.occ.addRectangle(x, y, 0, dr, dz)
        gmsh_ids.append(_id)

        y += dz

    if abs(Helix.z[1] - y) >= 0:
        _id = gmsh.model.occ.addRectangle(x, y, 0, dr, abs(Helix.z[1] - y))
        gmsh_ids.append(_id)

    # Now create air
    if AirData:
        r0_air = 0
        dr_air = Helix.r[1] * AirData[0]
        z0_air = Helix.z[0] * AirData[1]
        dz_air = abs(Helix.z[0] - Helix.z[1]) * AirData[1]
        _id = gmsh.model.occ.addRectangle(r0_air, z0_air, 0, dr_air, dz_air)

        ov, ovv = gmsh.model.occ.fragment([(2, _id)], [(2, i) for i in gmsh_ids])
        gmsh.model.occ.synchronize()
        return (gmsh_ids, (_id, dr_air, z0_air, dz_air))

    if debug:
        print(f"Helix/gmsh_ids: {gmsh_ids} ({len(gmsh_ids)})")

    # need to account for changes
    gmsh.model.occ.synchronize()
    return (gmsh_ids, ())


def gmsh_bcs(Helix: Helix, mname: str, ids: tuple, debug: bool = False) -> dict:
    """
    retreive ids for bcs in gmsh geometry
    """
    # print(f"Helix/gmsh_ids: ids={ids}")
    defs = {}
    (H_ids, Air_data) = ids

    prefix = ""
    if mname:
        prefix = f"{mname}_"

    # set physical name
    for i, id in enumerate(H_ids):
        ps = gmsh.model.addPhysicalGroup(2, [id])
        gmsh.model.setPhysicalName(2, ps, f"{prefix}Cu{i}")
        defs[f"{prefix}Cu{i}"] = ps

    # get BC ids
    gmsh.option.setNumber("Geometry.OCCBoundsUseStl", 1)

    zmin = Helix.z[0]
    zmax = Helix.z[1]

    bcs_defs = {
        f"{prefix}rInt": [Helix.r[0], zmin, Helix.r[0], zmax],
        f"{prefix}rExt": [Helix.r[1], zmin, Helix.r[1], zmax],
    }

    for key in bcs_defs:
        defs[key] = create_bcs(key, bcs_defs[key], 1)

    return defs
