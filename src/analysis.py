def get_closest_bead(
        atoms_df,
        prot,
        res,
        copy
):
    beads_df = atoms_df[(atoms_df["bead"] == 1) & (atoms_df["prot"] == prot) & (atoms_df["copy"] == copy)]
    closest_i = 0
    cur_res = beads_df.iloc[0, beads_df.columns.get_loc("res")]
    min_dist = abs(cur_res - res)
    for i in range(len(beads_df)):
        cur_res = beads_df.iloc[i, beads_df.columns.get_loc("res")]
        cur_dist = abs(cur_res - res)
        if cur_dist <= min_dist:
            closest_i = i
            min_dist = cur_dist

    return beads_df.iloc[closest_i, beads_df.columns.get_loc("res")]


def get_coords(
        atoms_df,
        prot,
        res,
        copy
):
    # Check if residue is contained:
    atoms_search = atoms_df[(atoms_df["prot"] == prot) & (atoms_df["res"] == res) & (atoms_df["copy"] == copy)]
    if len(atoms_search) > 0:
        return atoms_search.iloc[0, atoms_search.columns.get_loc("xyz")]
    else:
        beads_df = atoms_df[(atoms_df["bead"] == 1) & (atoms_df["prot"] == prot) & (atoms_df["copy"] == copy)]

        # If the model contains no 10-residue beads.
        if len(beads_df) == 0:
            return None

        closest_i = 0
        cur_res = beads_df.iloc[0, beads_df.columns.get_loc("res")]
        min_dist = abs(cur_res - res)
        for i in range(len(beads_df)):
            cur_res = beads_df.iloc[i, beads_df.columns.get_loc("res")]
            cur_dist = abs(cur_res - res)
            if cur_dist <= min_dist:
                closest_i = i
                min_dist = cur_dist

        if min_dist > 6:
            return None
        else:
            return beads_df.iloc[closest_i, atoms_search.columns.get_loc("xyz")]


def is_res_10(
        atoms_df,
        prot,
        res,
        copy
):
    atom = atoms_df[(atoms_df["prot"] == prot) & (atoms_df["res"] == res) & (atoms_df["copy"] == copy)]
    if len(atom) == 0:
        return True

    if atom.iloc[0, atom.columns.get_loc("bead")]:
        return True

    return False
