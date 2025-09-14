from pathlib import Path

import math


def z_on_sphere(x, y, radius):
    # Convert curvature from degrees per unit to radians per unit
    # curvature_radians = math.radians(curvature)
    # # Calculate the radius of the sphere
    # radius = 1 / (curvature_radians * math.tan(math.radians(1/2)))
    # print(radius)

    # Calculate the z value using the equation of the sphere
    z = math.sqrt(radius**2 - x**2 - y**2)

    return z


if __name__ == "__main__":
    # # Example usage
    # x = 0  # Example x coordinate
    # y = 0   # Example y coordinate
    # curvature = 1 / 200  # Curvature in degrees per unit

    # print(z_on_sphere(150, 150, 100000))
    # # print(z_on_sphere(200,200, curvature))

    membrane_file = Path(Path.home(), "Documents/mtorc2/data/pdb/etc/membrane_top_1.pdb")
    f = open(membrane_file, "w")

    # x_0, y_0, z_0 = 166.415, 166.428, 249.481

    # MTOR GLY 203
    x_0, y_0, z_0 = 166.9824, 166.9824, 87.550

    radius = 2.5
    atom_id = 1

    for x_offset in range(-45,45):
        for y_offset in range(-45,45):
            x = x_offset * radius*2 + x_0
            y = y_offset * radius*2 + y_0
            z = z_0
            pdb_entry = "ATOM{:>7}  CA  BEA A{:>4}     {:>7.3f} {:>7.3f} {:>7.3f}  1.00  5.99           C\n".format(atom_id, atom_id, x,y,z)
            f.write(pdb_entry)

            atom_id = atom_id + 1

    f.close()