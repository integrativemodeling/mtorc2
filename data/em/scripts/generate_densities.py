from pathlib import Path
import pandas as pd

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

    param_df = pd.read_csv(Path(Path.home(), "mtorc2/data/params/130.csv"))

    mols = dict()
    chain = "A"

    subunits = ["MSIN1", "AKT1"]
    for subunit in subunits:
        print(subunit)

        subunit_row_ids = list(param_df[param_df["subunit"] == subunit].index)
        color = param_df.iloc[subunit_row_ids[0], param_df.columns.get_loc("color")]
        chain = param_df.iloc[subunit_row_ids[0], param_df.columns.get_loc("model_chain")]

        mol = st.create_molecule(
            name=subunit,
            sequence=seqs[subunit],
            chain_id=chain
        )
        mols[subunit] = mol

        print(subunit_row_ids)
        for row_id in subunit_row_ids:
            name = param_df.iloc[row_id, param_df.columns.get_loc("name")]
            pdb_file = param_df.iloc[row_id, param_df.columns.get_loc("pdb_file")]
            pdb_file = Path(glob_data_dir, pdb_file)
            pdb_chain = param_df.iloc[row_id, param_df.columns.get_loc("pdb_chain")]
            start_pdb = param_df.iloc[row_id, param_df.columns.get_loc("start_pdb")]
            end_pdb = param_df.iloc[row_id, param_df.columns.get_loc("end_pdb")]
            offset = param_df.iloc[row_id, param_df.columns.get_loc("offset")]

            atom_res = mol.add_structure(
                pdb_fn=str(pdb_file),
                chain_id=pdb_chain,
                soft_check=True,
                res_range=(start_pdb, end_pdb),
                offset=offset,
                model_num=1
            )

            if name in ["CRIM", "KINASEN", "KINASEC"]:
                density_prefix = Path(comps_dir, name)
                mol.add_representation(
                    residues=atom_res,
                    density_residues_per_component=res_per_comp,
                    density_prefix=density_prefix,
                    density_force_compute=False,
                    resolutions=[1,10],
                    color=color
                )
            else:
                mol.add_representation(
                    residues=atom_res,
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