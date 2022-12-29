# ----------------------------------------------------------------------------
# <                               ProCare                                    >
# ----------------------------------------------------------------------------
# The MIT License (MIT)
#
# Copyright (c) 2020 Université de Strasbourg
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
# IN THE SOFTWARE.
# ----------------------------------------------------------------------------



"""
# apply a tranformation matrix no new coordinates from
# a previous ProCare alignment
# matrices are stored in ProCare output procare.tsv
"""



import numpy as np
from time import strftime, localtime
import os


def xtract_matrix(procare_file_, line_number_):
    with open(procare_file_, 'r') as f:
        data = f.read().split('\n')
    matrix_line = data[line_number_]
    matrix_cpnts = matrix_line.split('\t')[33:49]
    matrix_cpnts = [float(n) for n in matrix_cpnts]
    matrix = np.array([[n for n in matrix_cpnts[0:4]],
                      [n for n in matrix_cpnts[4:8]],
                      [n for n in matrix_cpnts[8:12]],
                      [n for n in matrix_cpnts[12:16]]])

    return matrix




def transform(mol2_, matrix_, ofile_):
    with open(mol2_, 'r') as f:
        mol2 = f.read()
        molecule_area = mol2.split('@<TRIPOS>MOLECULE\n')[1].split('@<TRIPOS>')[0]
        molecule_area = "@<TRIPOS>MOLECULE\n" + molecule_area
        atom_area = mol2.split('@<TRIPOS>ATOM\n')[1].split('@<TRIPOS>')[0]
        atoms = atom_area.split('\n')
        #print(atoms)
        del atoms[-1]
        new_atom_area = "@<TRIPOS>ATOM"
        for l in atoms:
            cols = l.split()
            coords = np.array([float(cols[2]), float(cols[3]), float(cols[4]), 1])
                # Umeyama Rt matrix [x, y, z, 1]
            matrix = np.array(matrix_)
            trans_coords = np.matmul(matrix, coords)
            #print(trans_coords, coords)
            new_atom_area += ("\n{:>7} {:<8} {:>9.4f} {:>9.4f} {:>9.4f} "
                              "{:<5} {:>5} {:<8} {:>9}".format(cols[0],
                                                     cols[1],
                                                     trans_coords[0],
                                                     trans_coords[1],
                                                     trans_coords[2],
                                                     cols[5],
                                                     cols[6],
                                                     cols[7],
                                                     cols[8]
                                                    ))
        new_atom_area += "\n"
        final_area = mol2.split('@<TRIPOS>BOND\n')[1]
        final_area = "@<TRIPOS>BOND\n" + final_area


    if ofile_[-5:] != '.mol2':
            ofile_ += '.mol2'
    basename = os.path.basename(ofile_)
    name = os.path.splitext(basename)[0]

    header = ""
    header += "# Modified by ProCare\n"
    header += "# Modification time: {}\n".format(
                            strftime("%a %d %b %Y %H:%M:%S", localtime()))
    header += "# Name: {}.mol2\n\n".format(name)

    with open(ofile_, 'w') as of:
        of.write('{}{}{}{}'.format(header, molecule_area, new_atom_area, final_area))

    print('output to {}'.format(ofile_))




if __name__ == '__main__':

    import argparse


    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fileprocare', type=str, required=True,
                help='procare matrix file, i.e. procare.tsv')
    parser.add_argument('-l', '--line', type=int, required=True,
                help='line number, eg. 1, where to extract matrix components. Starts with 0 = header')
    parser.add_argument('-a', '--appliedtomol2', type=str, required=True, nargs='+',
                help='mol2 file to apply transformation')
    parser.add_argument('--prefix', type=str, required=False,
                help='prefix to output files')
    args = parser.parse_args()


    matrix = xtract_matrix(args.fileprocare, args.line)
    for mol2 in args.appliedtomol2:
        filename = os.path.basename(mol2)
        prefix = 'rot'
        if args.prefix is not None:
            prefix = args.prefix
        ofile = '{}_{}'.format(prefix, filename)
        transform(mol2, matrix, ofile)
