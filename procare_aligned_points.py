
import numpy as np
from sklearn.neighbors import NearestNeighbors
import os
from time import strftime, localtime


def read_mol2(mol2_):
    with open(mol2_, 'r') as f:
        mol2 = f.read()

    atom_area = mol2.split('@<TRIPOS>ATOM\n')[1].split('@<TRIPOS>')[0]
    atoms_tmp = atom_area.split('\n')
    del atoms_tmp[-1]
    atoms = []
    for atm in atoms_tmp:
        if atm == '':
            continue
        atoms.append(atm.split())
    return atoms



def aligned_points(atoms_1_, atoms_2_):
    #print(len(atoms_1_))
    #print(len(atoms_2_))


    coords_1 = np.array([[float(p[2]), float(p[3]), float(p[4])] for p in atoms_1_])
    coords_2 = np.array([[float(p[2]), float(p[3]), float(p[4])] for p in atoms_2_])
    N = NearestNeighbors(n_neighbors=len(coords_1), radius=1.5, algorithm="ball_tree").fit(coords_1)
    distances, indices = N.radius_neighbors(np.array(coords_2))
    aligned_atoms_1 = []
    aligned_atoms_2 = []
    #print('len indices', len(indices))
    for i in range(len(coords_2)):
        #print('indices i', indices[i])
        for j in indices[i]:
            #print('j', j)
            if atoms_2_[i][1] == atoms_1_[j][1]:
                aligned_atoms_1.append(tuple(atoms_1_[j]))
                aligned_atoms_2.append(tuple(atoms_2_[i]))

    aligned_atoms_1 = list(set(aligned_atoms_1))
    aligned_atoms_2 = list(set(aligned_atoms_2))

    #print(len(aligned_atoms_1))
    #print(len(aligned_atoms_2))

    return aligned_atoms_1, aligned_atoms_2




def write_mol2(aligned_atoms_, ofile_):
    atm_to_resi = {'CA': 'GLY',
                    'CZ': 'PHE',
                    'O': 'ALA',
                    'OD1': 'ASP',
                    'OG': 'SER',
                    'N': 'ALA',
                    'NZ': 'LYS',
                    'DU': 'CUB',
                    }
    name = os.path.basename(ofile_)
    name = os.path.splitext(name)[0]
    of_string = ""
    of_string += "# Modified by ProCare\n"
    of_string += "# Modification time: {}\n".format(
                            strftime("%a %d %b %Y %H:%M:%S", localtime()))
    of_string += "# Name: {}.mol2\n\n".format(name)

    of_string += "@<TRIPOS>MOLECULE\n"
    of_string += "{}\n".format(name)
    of_string += "{:>5}{:>6}{:>6}{:>6}{:>6}\n".format(
                                            len(aligned_atoms_), 0, 0, 0, 0)
    of_string += "{}\n".format('PROTEIN')
    of_string += "NO_CHARGES\n"
    of_string += "@<TRIPOS>ATOM"

    for i, point in enumerate(aligned_atoms_):
        atm, x, y, z, atm_type, resn, resi, c =  point[1:]
        x = float(x)
        y = float(y)
        z = float(z)
        of_string += ("\n{:>7} {:<8} {:>9.4f} {:>9.4f} {:>9.4f} "
                      "{:<5} {:>5} {:<8} {:>9}".format(i+1, # replace original indice
                        atm, x, y, z, atm_type, i+1, atm_to_resi[atm]+str(i+1), c
                                                ))
    of_string += "\n@<TRIPOS>BOND"
    
    with open(ofile_, 'w') as of:
        of.write(of_string)
    print("written mol2 to {}".format(ofile_))




def ph4_contribution_score(aligned_atoms_1_, aligned_atoms_2_, atoms_1_, atoms_2_):
    if len(atoms_1_) < len(atoms_2_):
        ref_aligned_atoms = aligned_atoms_2_
        fit_aligned_atoms = aligned_atoms_1_
        ref_atoms = atoms_2_
        fit_atoms = atoms_1_
    else:
        ref_aligned_atoms = aligned_atoms_1_
        fit_aligned_atoms = aligned_atoms_2_
        ref_atoms = atoms_1_
        fit_atoms = atoms_2_

    ratio_CA_in_aligned = 0
    ratio_CZ_in_aligned = 0
    ratio_N_in_aligned = 0
    ratio_NZ_in_aligned = 0
    ratio_O_in_aligned = 0
    ratio_OD1_in_aligned = 0
    ratio_OG_in_aligned = 0
    ratio_DU_in_aligned = 0
    ratio_aligned = 0

    for _, atm, x, y, z, atm_type, resn, resi, c in fit_aligned_atoms:
        if atm == 'CA':
            ratio_CA_in_aligned += 1/len(fit_aligned_atoms)
        elif atm == 'CZ':
            ratio_CZ_in_aligned += 1/len(fit_aligned_atoms)
        elif atm == 'O':
            ratio_O_in_aligned += 1/len(fit_aligned_atoms)
        elif atm == 'OD1':
            ratio_OD1_in_aligned += 1/len(fit_aligned_atoms)
        elif atm == 'OG':
            ratio_OG_in_aligned += 1/len(fit_aligned_atoms)
        elif atm == 'N':
            ratio_N_in_aligned += 1/len(fit_aligned_atoms)
        elif atm == 'NZ':
            ratio_NZ_in_aligned += 1/len(fit_aligned_atoms)
        elif atm == 'DU':
            ratio_DU_in_aligned += 1/len(fit_aligned_atoms)
        ratio_aligned += 1/len(fit_atoms)

    ratio_aligned = round(ratio_aligned, 4)
    ratio_CA_in_aligned = round(ratio_CA_in_aligned, 4)
    ratio_CZ_in_aligned = round(ratio_CZ_in_aligned, 4)
    ratio_O_in_aligned = round(ratio_O_in_aligned, 4)
    ratio_N_in_aligned = round(ratio_N_in_aligned, 4)
    ratio_OD1_in_aligned = round(ratio_OD1_in_aligned, 4)
    ratio_OG_in_aligned = round(ratio_OG_in_aligned, 4)
    ratio_NZ_in_aligned = round(ratio_NZ_in_aligned, 4)
    ratio_DU_in_aligned = round(ratio_DU_in_aligned, 4)


    return ratio_aligned, \
            ratio_CA_in_aligned, ratio_CZ_in_aligned, \
            ratio_N_in_aligned, ratio_NZ_in_aligned, \
            ratio_O_in_aligned, ratio_OD1_in_aligned, \
            ratio_OG_in_aligned, ratio_DU_in_aligned




if __name__ == '__main__':
    
    import argparse

    """ output aligned points in two IChem VolSite cavities
        ease visualization """



    parser = argparse.ArgumentParser()
    parser.add_argument('-c1', '--cav1', type=str, required=True,
                help='cavity1')
    parser.add_argument('-c2', '--cav2', type=str, required=True,
                help='cavity2')
    parser.add_argument('-o1', '--ocav1', type=str, required=True,
                help='output aligned cavity1')
    parser.add_argument('-o2', '--ocav2', type=str, required=True,
                help='output aligned cavity2')
    parser.add_argument('--contrib', action='store_true', required=False,
                help='output individual proportions of aligned ph4')
    args = parser.parse_args()


    atoms_1 = read_mol2(args.cav1)
    atoms_2 = read_mol2(args.cav2)

    aligned_atoms_1, aligned_atoms_2 = aligned_points(atoms_1, atoms_2)
    ratio_aligned, ratio_CA_in_aligned, ratio_CZ_in_aligned, \
    ratio_N_in_aligned, ratio_NZ_in_aligned, ratio_O_in_aligned, \
    ratio_OD1_in_aligned, ratio_OG_in_aligned, \
        ratio_DU_in_aligned = ph4_contribution_score(
            aligned_atoms_1, aligned_atoms_2, atoms_1, atoms_2)

    # output aligned points
    write_mol2(aligned_atoms_1, args.ocav1)
    write_mol2(aligned_atoms_2, args.ocav2)


    if args.contrib is not None:

        if not os.path.isfile('procare_scores_contribution.tsv'):
            with open('procare_scores_contribution.tsv', 'w') as of:
                of.write('cavity1\tcavity2\tratio_aligned\tCA_contrib\tCZ_contrib\t'
                        'O_contrib\tN_contrib\tOD1_contrib\tOG_contrib\t'
                        'NZ_contrib\tDU_contrib\n')

        with open('procare_scores_contribution.tsv', 'a') as of:
            of.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                os.path.basename(args.cav1),
                os.path.basename(args.cav2),
                ratio_aligned,
                ratio_CA_in_aligned,
                ratio_CZ_in_aligned,
                ratio_O_in_aligned,
                ratio_N_in_aligned,
                ratio_OD1_in_aligned,
                ratio_OG_in_aligned,
                ratio_NZ_in_aligned,
                ratio_DU_in_aligned))


