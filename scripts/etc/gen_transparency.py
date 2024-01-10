from pathlib import Path
import pandas as pd


def gen_transparency(
        pdb_file,
        out_file,
        coord,
        view,
        size,
        omit_chains
):
    f = open(pdb_file, "r")

    atom_dict = {"chain": list(), "res_id": list(), "res_type": list(), "x": list(), "y": list(), "z": list()}

    bead_dict = {"chain": list(), "res_id": list(), "x": list(), "y": list(), "z": list()}

    for x in f:
        if x[0:4] == "ATOM":
            res_id = int(x[22:26])
            coords = float(x[31:38]), float(x[39:46]), float(x[47:54])
            chain = x[21]
            res_type = x[17:20]
            if chain in omit_chains:
                continue

            if res_type == "BEA":
                bead_dict["chain"].append(chain)
                bead_dict["res_id"].append(res_id)
                bead_dict["x"].append(coords[0])
                bead_dict["y"].append(coords[1])
                bead_dict["z"].append(coords[2])

            atom_dict["chain"].append(chain)
            atom_dict["res_id"].append(res_id)
            atom_dict["res_type"].append(res_type)
            atom_dict["x"].append(coords[0])
            atom_dict["y"].append(coords[1])
            atom_dict["z"].append(coords[2])

    atom_df = pd.DataFrame(atom_dict)
    # atom_df["x"] = atom_df["x"] - atom_df["x"].min()
    # atom_df["y"] = atom_df["y"] - atom_df["y"].min()
    # atom_df["z"] = atom_df["z"] - atom_df["z"].min()

    avg_coord_dict = {"chain": list(), "start": list(), "end": list(), "avg": list()}
    f = open(out_file, "w")
    chains = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"]
    for omit_chain in omit_chains:
        chains.remove(omit_chain)
    for chain_id in chains:
        chain_df = atom_df[atom_df["chain"] == chain_id]

        for i in range(chain_df["res_id"].max() // size + 1):
            section_df = chain_df[(chain_df.res_id >= (i*size)) & (chain_df.res_id < ((i+1)*size))]

            # break

            if len(section_df) == 0:
                continue

            avg = 0
            for j in range(len(section_df)):
                coord_val = section_df.iloc[j, section_df.columns.get_loc(coord)]

                avg = avg + coord_val

            avg = avg / len(section_df)
            start = section_df["res_id"].min()
            end = section_df["res_id"].max()

            avg_coord_dict["chain"].append(chain_id)
            avg_coord_dict["start"].append(start)
            avg_coord_dict["end"].append(end)
            avg_coord_dict["avg"].append(avg)

    avg_coord_df = pd.DataFrame(avg_coord_dict)

    avg_coord_df["avg"] = avg_coord_df["avg"] - avg_coord_df["avg"].min()
    max = avg_coord_df["avg"].max()
    min = avg_coord_df["avg"].min()
    print(max, min)

    for i in range(len(avg_coord_df)):
        chain_id, start, end, avg_coord = avg_coord_df.iloc[i,0:4]
        print(chain_id, start, end, avg_coord)
        # trans = ((avg_coord - min) / (max - min) * 21.55)**1.5

        if view == 0:
            trans = ((max - avg_coord) / (max - min))**(1/2) * 100
        else:
            trans = ((avg_coord - min) / (max - min))**(1/2) * 100

        # trans = ((avg_coord - min) / (max - min)) * 100
        # trans = ((avg_coord - min) / (max - min))**2 * 100
        # trans = 100 - (avg_z - min)**2 / (max - min)
        # trans = 100 - int(trans)
        # else:
        #     trans = int(trans)

        f.write("transparency #1/{}:{}-{} {} ribbons\n".format(chain_id, start, end, trans))

    bead_df = pd.DataFrame(bead_dict)
    bead_df["x"] = bead_df["x"] - atom_df["x"].min()
    bead_df["y"] = bead_df["y"] - atom_df["y"].min()
    bead_df["z"] = bead_df["z"] - atom_df["z"].min()

    print(min, max)
    for i in range(len(bead_df)):
        chain_id, res_id, x, y, z = bead_df.iloc[i,0:5]
        coords = dict()
        coords["x"], coords["y"], coords["z"] = x, y, z

        if coords[coord] > max:
            coords[coord] = max
        elif coords[coord] < min:
            coords[coord] = min

        if view == 0:
            trans = ((max - coords[coord]) / (max - min))**(1/2) * 100
        else:
            trans = ((coords[coord] - min) / (max - min))**(1/2) * 100

        if res_id == 1175:
            print(coords, trans)

        f.write("transparency #1/{}:{} {} atoms\n".format(chain_id, res_id, trans))

    f.close()


if __name__ == "__main__":
    for view in [0,1]:
        for coord in ["x", "y", "z"]:
            gen_transparency(
                pdb_file=Path(Path.home(), "Documents/mtorc2/data/models/submission_3//136_2_1/cluster.0/cluster_center_model.pdb"),
                out_file=Path(Path.home(), "Documents/mtorc2/scripts/chimerax/transparency_{}_{}.cxc".format(coord, view)),
                coord=coord,
                view=view,
                size=30,
                omit_chains=["E", "J"]
            )
