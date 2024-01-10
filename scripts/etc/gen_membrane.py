from pathlib import Path

if __name__ == "__main__":
    membrane_file = Path("./membrane_1.pdb")
    f = open(membrane_file, "w")

    x_0, y_0, z_0 = 166.415, 166.428, 237.481
    # x_0, y_0, z_0 = 166.415, 166.428, 242.481

    radius = 2.5
    atom_id = 1

    for x_offset in range(-40,40):
        for y_offset in range(-40,40):
            x = x_offset * radius*2 + x_0
            y = y_offset * radius*2 + y_0
            z = z_0
            pdb_entry = "ATOM{:>7}  CA  BEA A{:>4}     {:>7.3f} {:>7.3f} {:>7.3f}  1.00  5.99           C\n".format(atom_id, atom_id, x,y,z)
            f.write(pdb_entry)

            atom_id = atom_id + 1

    f.close()