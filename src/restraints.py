import IMP
import IMP.pmi
import IMP.pmi.io.crosslink
import IMP.pmi.restraints
import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.stereochemistry
import ihm.cross_linkers
from pathlib import Path


def get_ev_restraints(
        mols,
        clones=None
):
    rs = list()

    r_1 = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
        included_objects=[mols[component] for component in mols.keys()],
        resolution=10
    )
    rs.append(r_1)

    components = list(mols.keys())
    if clones:
        for component in mols.keys():
            for component_2 in components:
                r_2 = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
                    included_objects=[mols[component], clones[component_2]],
                    resolution=10
                )
                rs.append(r_2)
            components.remove(component)

    return rs


def get_crs(
   mols,
):
    crs = []
    for component in mols.keys():
        mol = mols[component]
        cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mol)
        crs.append(cr)

    return crs


def get_xl_restraint(
        hier,
        xl_file,
        w,
        type
):
    xl_dir = Path(Path.home(), "mtorc2/data/xlms")
    xldbkwc = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
    xldbkwc.set_protein1_key("prot1")
    xldbkwc.set_protein2_key("prot2")
    xldbkwc.set_residue1_key("res1")
    xldbkwc.set_residue2_key("res2")

    xldb = IMP.pmi.io.crosslink.CrossLinkDataBase(xldbkwc)
    xldb.create_set_from_file(str(Path(xl_dir, xl_file)))

    if type == "DSS":
        linker=ihm.cross_linkers.dss
    elif type == "EDC":
        linker=ihm.cross_linkers.edc
    else:
        raise RuntimeError("Specify XL type")

    r = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
        root_hier=hier,
        database=xldb,
        length=35.,
        slope=0.01,
        resolution=1.,
        label=type,
        linker=linker,
        weight=w
    )

    return r