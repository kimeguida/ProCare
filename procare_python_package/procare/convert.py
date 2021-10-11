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


"""Conversion between mol2 and pcd format of protein IChem VolSite cavities"""


import os
from time import strftime, localtime
import numpy as np


class _mol2_:

    def _mol2_to_pcd(self, ifile_, color_):
        """ Extracts coordinates from mol2 files and convert into pcd format """

        if ifile_[-5:] != '.mol2':
            print("incorrect file extension")
            print("file format may be wrong --> no output")

        try:
            with open(ifile_, "r") as f:
                mol2 = f.read()

            atom_area = mol2.split('@<TRIPOS>ATOM\n')[1].split('@<TRIPOS>')[0]
            atoms_tmp = atom_area.split('\n')
            #print(atoms)
            del atoms_tmp[-1]
            atoms = []
            for atm in atoms_tmp:
                if atm == '':
                    continue
                atoms.append(atm)

        except:
            print("Cannot process mol2 {}".format(ifile_))
            return -1, None, None


        ofilename = os.path.basename(ifile_).replace("mol2", "pcd")    
        try:
            ofile = open(ofilename, "w")
        except IOError:
            print("mol2_to_pcd: Cannot write to current directory: "
                  "{}. Please check for access rights.".format(os.getcwd()))
            ofile.close()
            return -1, None, None

        properties = []
        colors = []
        ofile.write("VERSION .7\nFIELDS x y z rgb\nSIZE 4 4 4 4\n"
                    "TYPE F F F F\nCOUNT 1 1 1 1\n")
        ofile.write("WIDTH {0}\nHEIGHT 1\nVIEWPOINT 0 0 0 1 0 0 0\n"
                    "POINTS {0}\nDATA ascii".format(len(atoms)))
        #print(atoms)
        for atm in atoms:
            cols = atm.split()
            index = int(cols[0])
            atom = str(cols[1])
            x = float(cols[2])
            y = float(cols[3])
            z = float(cols[4])
            properties.append([index, atom])
            colors.append(color_[atom])
            ofile.write("\n{} {} {} {}".format(x, y, z, color_[atom]))

        ofile.close()                                                    
        
        #print(ofilename)
        return ofilename, properties, colors



class _pcd_:

    def _get_coordinates(self, ifile_):

        if ifile_[-4:] != '.pcd':
            print("incorrect file extension")
            print("file format may be wrong --> no output")
            if "DATA ascii" not in open(ifile_, 'r').read():
                return -1

        coordinates = []
        errors = []
        with open(ifile_, 'r') as f:
            data = f.read().split('\n')
        data = [l for l in data[10:] if l != '']

        for line_ in data:
            x, y, z, rgb = line_.split()
            try: 
                x = float(x)
            except ValueError:
                errors.append(-1)
            try:
                y = float(y)
            except ValueError:
                errors.append(-1)
            try:
                z = float(z)
            except ValueError:
                errors.append(-1)
            try:
                rgb = int(rgb)
            except ValueError:
                errors.append(-1)

            if -1 in errors:
                print("Errors while parsing PCD")
                break

            coordinates.append([x, y, z, rgb])

        if -1 not in errors:
            self.coordinates = coordinates
            return coordinates
        else:
            return -1



    def _write_mol2(self, ofile_, coordinates_, atom_, atom_type_, residue_, 
                                                    macromol_="PROTEIN"):

        
        if ofile_[-5:] != '.mol2':
            ofile_ += '.mol2'
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
                                                len(coordinates_), 0, 0, 0, 0)
        of_string += "{}\n".format(macromol_)
        of_string += "NO_CHARGES\n"
        of_string += "@<TRIPOS>ATOM"

        for i, point in enumerate(coordinates_):
            x, y, z, rgb = [*point]
            of_string += ("\n{:>7} {:<8} {:>9.4f} {:>9.4f} {:>9.4f} "
                          "{:<5} {:>5} {:<8} {:>9}".format(i+1,
                                                     atom_[rgb],
                                                     x,
                                                     y,
                                                     z,
                                                     atom_type_[rgb],
                                                     i+1,
                                                     residue_[rgb]+str(i+1),
                                                     0.0000
                                                    ))
        of_string += "\n@<TRIPOS>BOND"
        
        with open(ofile_, 'w') as of:
            of.write(of_string)
        print("written mol2 to {}".format(ofile_))
        self.type = "pcd"
        return ofile_

    
    def _pcd_to_mol2(self, ifile_, atom_, atom_type_, residue_, 
                                                    macromol_="PROTEIN"):
                
        coordinates = self._get_coordinates(ifile_)
        if coordinates != -1:
            ofile = os.path.basename(ifile_).replace('pcd', 'mol2')
            if self._write_mol2(ofile, coordinates, atom_, atom_type_, 
                                                residue_, macromol_) == ofile:
                self.ifile = ifile_
                return ofile
            else:
                return -1


    

class _volsite_cavity_(_mol2_, _pcd_):
    
    def __init__(self):

        __COLOR = {"OG":8204959,
                   "N":30894,
                   "O":15219528,
                   "NZ":15231913,
                   "CZ":4646984,
                   "CA":16741671,
                   "DU":7566712,
                   "OD1":0,}


        __ATOM = {val:key for key, val in __COLOR.items()}


        __ATOM_TYPE = {"OG":"O.3",
                          "N":"N.am",
                          "O":"O.2",
                          "NZ":"N.4",
                          "CZ":"C.ar",
                          "CA":"C.3",
                          "DU":"H",
                          "OD1":"O.co2",}

        __RESIDUE = {"OG":"SER",
                        "N":"ALA",
                        "O":"ALA",
                        "NZ":"LYS",
                        "CZ":"PHE",
                        "CA":"GLY",
                        "DU":"CUB",
                        "OD1":"ASP",}

        self.COLOR = __COLOR

        self.ATOM = __ATOM

        self.ATOM_TYPE = {key:__ATOM_TYPE[val] 
                                for key, val in __ATOM.items()}
        self.RESIDUE = {key:__RESIDUE[val] 
                                for key, val in __ATOM.items()}

        

    def mol2_to_pcd(self, ifile_):
        return self._mol2_to_pcd(ifile_, self.COLOR)
        



    def pcd_to_mol2(self, ifile_):
        return self._pcd_to_mol2(ifile_, self.ATOM,
                                         self.ATOM_TYPE,
                                         self.RESIDUE)



    def write_mol2(self, ofile_, coordinates_):
        return self._write_mol2(ofile_, coordinates_, self.ATOM, self.ATOM_TYPE,
                    self.RESIDUE)



if __name__ == '__main__':
    
    import argparse
    import os
    from time import strftime, localtime
    import numpy as np

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, 
                help="input file",
                required=True)
    parser.add_argument('-t', '--itype', type=str, 
                help="input file type: mol2 or pcd",
                choices=["mol2", "pcd"],
                required=True)
    parser.add_argument('-m', '--macromol', type=str, 
                help="macromolecule type: cavity, ...",
                choices=["cav"],
                required=True,
                default="cav")

    args = parser.parse_args()

    if args.macromol == "cav":
        molecule = _volsite_cavity_()
        if args.itype == "pcd":
            molecule.pcd_to_mol2(args.input)
        elif args.itype == "mol2":
            molecule.mol2_to_pcd(args.input)

        #coords = [[1, 2, 3, 8204959], [4, 5, 6, 8204959]]
        #molecule.write_mol2('test.mol2', coords)