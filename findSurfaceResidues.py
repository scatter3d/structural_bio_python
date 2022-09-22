    #    Finds those atoms on the surface of a protein
    #    that have at least 'cutoff' exposed A**2 surface area.
    #    import script to pymol using GUI
    #    USAGE - findSurfaceAtoms [ selection, [ cutoff, [ doShow ]]]
    #    ex. findSurfaceResidues all, 2.5, doshow=0(optional)


# import 
from __future__ import print_function
from pymol import cmd


def findSurfaceAtoms(selection="all", cutoff=2.5, quiet=1):

    cutoff, quiet = float(cutoff), int(quiet)

    tmpObj = cmd.get_unused_name("_tmp")
    cmd.create(tmpObj, "(" + selection + ") and polymer", zoom=0)

    cmd.set("dot_solvent", 1, tmpObj)
    cmd.get_area(selection=tmpObj, load_b=1)

    # threshold on what one considers an "exposed" atom (in A**2):
    cmd.remove(tmpObj + " and b < " + str(cutoff))

    selName = cmd.get_unused_name("exposed_atm_")
    cmd.select(selName, "(" + selection + ") in " + tmpObj)

    cmd.delete(tmpObj)

    if not quiet:
        print("Exposed atoms are selected in: " + selName)

    return selName


def findSurfaceResidues(selection="all", cutoff=2.5, doShow=0, quiet=1):

    cutoff, doShow, quiet = float(cutoff), int(doShow), int(quiet)

    selName = findSurfaceAtoms(selection, cutoff, quiet)

    exposed = set()
    cmd.iterate(selName, "exposed.add((chain,resv))", space=locals())

    selNameRes = cmd.get_unused_name("exposed_res_")
    cmd.select(selNameRes, "byres " + selName)

    if not quiet:
        print("Exposed residues are selected in: " + selNameRes)

    if doShow:
        cmd.show_as("spheres", "(" + selection + ") and polymer")
        cmd.color("white", selection)
        cmd.color("yellow", selNameRes)
        cmd.color("red", selName)

    return sorted(exposed)

cmd.extend("findSurfaceAtoms", findSurfaceAtoms)
cmd.extend("findSurfaceResidues", findSurfaceResidues)
