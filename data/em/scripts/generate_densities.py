from pathlib import Path
import IMP
import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros


def get_densities(
        comps_dir,
        res_per_comp
):
    m = IMP.Model()
    s = IMP.pmi.topology.System(m)

    glob_data_dir = Path(Path.home(), "mtorc2/data")
    fasta_file = Path(glob_data_dir, "fasta/mtorc2.fasta")
    seqs = IMP.pmi.topology.Sequences(fasta_fn=str(fasta_file))
    st = s.create_state()

    subunits = ["MSIN1", "AKT1"]
    structures = dict()
    pdb_dir = Path(glob_data_dir, "em/models/J252")
    structures["MSIN1"] = list()
    structures["MSIN1"].append((Path(pdb_dir, "mSIN1_CRIM_AF-Q9BPZ7-F1_conf02.pdb"), "A", "red", 158, 267, 158, 267, 0, "CRIM"))
    structures["AKT1"] = list()
    structures["AKT1"].append((Path(pdb_dir, "AKT1_N-lobe_AF-P31749-F1_conf02.pdb"), "A", "purple", 146, 229, 146, 229, 0, "KINASEN"))
    structures["AKT1"].append((Path(pdb_dir, "AKT1_C-lobe_AF-P31749-F1_conf02.pdb"), "A", "purple", 230, 480, 234, 425, 0, "KINASEC"))

    mols = dict()
    chain = "A"
    for subunit in subunits:
        print(subunit)

        mol = st.create_molecule(
            name=subunit,
            sequence=seqs[subunit],
            chain_id=chain
        )
        mols[subunit] = mol

        for file, struct_chain, color, start, end, start_struct, end_struct, offset, prefix in structures[subunit]:
            print(file, struct_chain, color)
            atom_res = mol.add_structure(
                pdb_fn=str(file),
                chain_id=struct_chain,
                soft_check=True,
                res_range=(start_struct, end_struct),
                offset=offset,
                model_num=1
            )

            if prefix:
                density_prefix = Path(comps_dir, prefix)
                mol.add_representation(
                    residues=atom_res,
                    density_residues_per_component=res_per_comp,
                    density_prefix=str(density_prefix),
                    density_force_compute=True,
                    resolutions=[1,res_per_comp,10],
                    color=color
                )

            else:
                mol.add_representation(
                    residues=mol.get_atomic_residues(),
                    resolutions=[1,10],
                    color=color
                )

    s.build()



if __name__ == "__main__":
    comps_dir = Path(Path.home(), "mtorc2/data/em/comps//130")
    get_densities(
        comps_dir=comps_dir,
        res_per_comp=10
    )