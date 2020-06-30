# ----------------------------------------------------------------------------
# <                               ProCare                                    >
# ----------------------------------------------------------------------------
# The MIT License (MIT)
#
# Copyright (c) 2020 Merveille Eguida
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




class _similarity_metrics_:

    

    def tanimoto_similarity(self): # mettre en arg la classe score et utuliser les instances ...   
        similarity = float(self.n_identity)/\
            (self.fitsize + self.refsize - self.n_identity)
        return round(similarity, 4)

    
    def tversky_similarity(self, alpha_=0.95, beta_=0.05):
        similarity = float(self.n_identity)/(alpha_ * (self.fitsize - self.n_identity)\
                                        + beta_ * (self.refsize - self.n_identity)\
                                        + self.n_identity)
        return round(similarity, 4)

    
    def cosine_similarity(self):
        import math
        similarity = float(self.n_identity)/(math.sqrt(self.fitsize)*math.sqrt(self.refsize))
        return round(similarity, 4)

    
    def dice_similarity(self):
        similarity = 2 * float(self.n_identity)/(self.fitsize + self.refsize)
        return round(similarity, 4)


    def per_score(self):
        return round(float(self.n_identity)/self.fitsize, 4)

    
    def weighted_score(self):

        frequency = {"CA": 0.3817,
                     "CZ": 0.1110,
                     "O": 0.0659,
                     "OG": 0.0867,
                     "OD1": 0.0593,
                     "N": 0.1376,
                     "NZ": 0.0579,
                     "DU": 0.0999}

        similarity = float(self.CA)/self.fitsize * 1/frequency["CA"] +\
                         float(self.CZ)/self.fitsize * 1/frequency["CZ"] +\
                         float(self.N)/self.fitsize * 1/frequency["N"] +\
                         float(self.OD1)/self.fitsize * 1/frequency["OD1"] +\
                         float(self.NZ)/self.fitsize * 1/frequency["NZ"] +\
                         float(self.OG)/self.fitsize * 1/frequency["OG"] +\
                         float(self.O)/self.fitsize * 1/frequency["O"] +\
                         float(self.DU)/self.fitsize * 1/frequency["DU"]
        return round(similarity, 4)





class _distances_:

    def hamming_distance(self):
        distance = self.fitsize + self.refsize - (2 * self.n_identity)
        return round(distance, 4)


    def soergel_distance(self):
        distance = float(self.fitsize + self.refsize- (2 * self.n_identity))/\
                                (self.fitsize + self.refsize - self.n_identity)
        return round(distance, 4)



    

class _ph4_strict_(_similarity_metrics_, _distances_):
    """ 1-NN search 
        and distance < D 
        and strict correspondence of properties """
    


    def __init__(self, source_coordinates_, target_coordinates_, 
                 source_properties_, target_properties_, distance_threshold_):

        import numpy as np
        from sklearn.neighbors import NearestNeighbors
        
        self.source_coordinates = source_coordinates_
        self.target_coordinates = target_coordinates_
        self.source_properties = source_properties_
        self.target_properties = target_properties_
        self.distance_threshold = distance_threshold_
        self.n_identity = 0
        self.n_different = 0
        self.CA = 0
        self.CZ = 0
        self.N = 0
        self.NZ = 0
        self.O = 0
        self.OG = 0
        self.OD1 = 0
        self.DU = 0

        if len(self.source_coordinates) > len(self.target_coordinates):
            N = NearestNeighbors(n_neighbors=1, 
                                 algorithm="ball_tree").fit(self.source_coordinates)
            self.distances, self.indices = N.kneighbors(np.array(self.target_coordinates))
            self.fitProp = self.target_properties
            self.refProp = self.source_properties
        
        else:
            N = NearestNeighbors(n_neighbors=1, 
                                 algorithm="ball_tree").fit(self.target_coordinates)
            self.distances, self.indices = N.kneighbors(np.array(self.source_coordinates))
            self.fitProp = self.source_properties
            self.refProp = self.target_properties

        self.fitsize = len(self.fitProp)
        self.refsize = len(self.refProp)


    def get_similarity_by_rules(self):
        
        for i in range(len(self.fitProp)):
            if self.distances[i][0] <= self.distance_threshold:
                if self.fitProp[i][0] == i+1: # if indice corresponds to indice in mol2 file ## check

                    if (self.fitProp[i][1] == "CA") and\
                            (self.refProp[self.indices[i][0]][1] == "CA"):
                        self.n_identity += 1
                        self.CA += 1

                    elif (self.fitProp[i][1] == "CZ") and\
                            (self.refProp[self.indices[i][0]][1] == "CZ"):
                        self.n_identity += 1
                        self.CZ += 1
 
                    elif (self.fitProp[i][1] == "N") and\
                            (self.refProp[self.indices[i][0]][1] == "N"):
                        self.n_identity += 1
                        self.N += 1
                        
                    elif (self.fitProp[i][1] == "O") and\
                            (self.refProp[self.indices[i][0]][1] == "O"):
                        self.n_identity += 1
                        self.O += 1

                    elif (self.fitProp[i][1] == "NZ") and\
                            (self.refProp[self.indices[i][0]][1] == "NZ"):
                        self.n_identity += 1
                        self.NZ += 1

                    elif (self.fitProp[i][1] == "OD1") and\
                            (self.refProp[self.indices[i][0]][1] == "OD1"):
                        self.n_identity += 1
                        self.OD1 += 1

                    elif (self.fitProp[i][1] == "OG") and\
                            (self.refProp[self.indices[i][0]][1] == "OG"):
                        self.n_identity += 1
                        self.OG += 1

                    elif (self.fitProp[i][1] == "DU") and\
                            (self.refProp[self.indices[i][0]][1] == "DU"):
                        self.n_identity += 1
                        self.DU += 1

        ratio_aligned = round(float(self.n_identity)/self.fitsize, 4)
        
        if self.n_identity != 0:
            ratio_CA_in_aligned = round(float(self.CA)/self.n_identity, 4)
            ratio_CZ_in_aligned = round(float(self.CZ)/self.n_identity, 4)
            ratio_N_in_aligned = round(float(self.N)/self.n_identity, 4)
            ratio_NZ_in_aligned = round(float(self.NZ)/self.n_identity, 4)
            ratio_O_in_aligned = round(float(self.O)/self.n_identity, 4)
            ratio_OD1_in_aligned = round(float(self.OD1)/self.n_identity, 4)
            ratio_OG_in_aligned = round(float(self.OG)/self.n_identity, 4)
            ratio_DU_in_aligned = round(float(self.DU)/self.n_identity, 4)
        else:
            ratio_CA_in_aligned = 0
            ratio_CZ_in_aligned = 0
            ratio_N_in_aligned = 0
            ratio_NZ_in_aligned = 0
            ratio_O_in_aligned = 0
            ratio_OD1_in_aligned = 0
            ratio_OG_in_aligned = 0
            ratio_DU_in_aligned = 0

        return ratio_aligned, \
            ratio_CA_in_aligned, ratio_CZ_in_aligned, \
            ratio_N_in_aligned, ratio_NZ_in_aligned, \
            ratio_O_in_aligned, ratio_OD1_in_aligned, \
            ratio_OG_in_aligned, ratio_DU_in_aligned

    



class _ph4_rules_(_similarity_metrics_, _distances_):
    """ 1-NN search 
        and distance < D 
        and correspondence of properties according to compatibility rules """

    def __init__(self, source_coordinates_, target_coordinates_, 
                 source_properties_, target_properties_, distance_threshold_):

        import numpy as np
        from sklearn.neighbors import NearestNeighbors
        
        self.source_coordinates = source_coordinates_
        self.target_coordinates = target_coordinates_
        self.source_properties = source_properties_
        self.target_properties = target_properties_
        self.distance_threshold = distance_threshold_
        self.n_identity = 0
        self.n_different = 0
        self.CA = 0
        self.CZ = 0
        self.N = 0
        self.NZ = 0
        self.O = 0
        self.OG = 0
        self.OD1 = 0
        self.DU = 0

        if len(self.source_coordinates) > len(self.target_coordinates):
            N = NearestNeighbors(n_neighbors=1, 
                                 algorithm="ball_tree").fit(self.source_coordinates)
            self.distances, self.indices = N.kneighbors(np.array(self.target_coordinates))
            self.fitProp = self.target_properties
            self.refProp = self.source_properties
        
        else:
            N = NearestNeighbors(n_neighbors=1, 
                                 algorithm="ball_tree").fit(self.target_coordinates)
            self.distances, self.indices = N.kneighbors(np.array(self.source_coordinates))
            self.fitProp = self.source_properties
            self.refProp = self.target_properties

        self.fitsize = len(self.fitProp)
        self.refsize = len(self.refProp)


    def get_similarity_by_rules(self):
        
        for i in range(len(self.fitProp)):
            if self.distances[i][0] <= self.distance_threshold:
                
                if (self.fitProp[i][1] == "CA") and\
                        (self.refProp[self.indices[i][0]][1] == "CA" or "CZ"):
                    self.n_identity += 1
                    self.CA += 1

                elif (self.fitProp[i][1] == "CZ") and\
                        (self.refProp[self.indices[i][0]][1] == "CA" or "CZ"):
                    self.n_identity += 1
                    self.CZ += 1

                elif (self.fitProp[i][1] == "N") and\
                        (self.refProp[self.indices[i][0]][1] == "N" or "NZ" or "OG"):
                    self.n_identity += 1
                    self.N += 1
                    
                elif (self.fitProp[i][1] == "O") and\
                        (self.refProp[self.indices[i][0]][1] == "O" or "OD1" or "OG"):
                    self.n_identity += 1
                    self.O += 1

                elif (self.fitProp[i][1] == "NZ") and\
                        (self.refProp[self.indices[i][0]][1] == "NZ" or "N" or "OG"):
                    self.n_identity += 1
                    self.NZ += 1

                elif (self.fitProp[i][1] == "OD1") and\
                        (self.refProp[self.indices[i][0]][1] == "OD1" or "O" or "OG"):
                    self.n_identity += 1
                    self.OD1 += 1

                elif (self.fitProp[i][1] == "OG") and\
                        (self.refProp[self.indices[i][0]][1] == "OG" or "N" or "O" or "NZ" or "OD1"):
                    self.n_identity += 1
                    self.OG += 1

                elif (self.fitProp[i][1] == "DU") and\
                        (self.refProp[self.indices[i][0]][1] == "DU"):
                    self.n_identity += 1
                    self.DU += 1

        ratio_aligned = round(float(self.n_identity)/self.fitsize, 4)
        
        if self.n_identity != 0:
            ratio_CA_in_aligned = round(float(self.CA)/self.n_identity, 4)
            ratio_CZ_in_aligned = round(float(self.CZ)/self.n_identity, 4)
            ratio_N_in_aligned = round(float(self.N)/self.n_identity, 4)
            ratio_NZ_in_aligned = round(float(self.NZ)/self.n_identity, 4)
            ratio_O_in_aligned = round(float(self.O)/self.n_identity, 4)
            ratio_OD1_in_aligned = round(float(self.OD1)/self.n_identity, 4)
            ratio_OG_in_aligned = round(float(self.OG)/self.n_identity, 4)
            ratio_DU_in_aligned = round(float(self.DU)/self.n_identity, 4)
        else:
            ratio_CA_in_aligned = 0
            ratio_CZ_in_aligned = 0
            ratio_N_in_aligned = 0
            ratio_NZ_in_aligned = 0
            ratio_O_in_aligned = 0
            ratio_OD1_in_aligned = 0
            ratio_OG_in_aligned = 0
            ratio_DU_in_aligned = 0

        return ratio_aligned, \
            ratio_CA_in_aligned, ratio_CZ_in_aligned, \
            ratio_N_in_aligned, ratio_NZ_in_aligned, \
            ratio_O_in_aligned, ratio_OD1_in_aligned, \
            ratio_OG_in_aligned, ratio_DU_in_aligned





class _ph4_ext_(_similarity_metrics_, _distances_):
    """  neighbors with distance < D 
         and strict correspondence of properties """

    
    def __init__(self, source_coordinates_, target_coordinates_, 
                 source_properties_, target_properties_, distance_threshold_):

        import numpy as np
        from sklearn.neighbors import NearestNeighbors

        self.source_coordinates = source_coordinates_
        self.target_coordinates = target_coordinates_
        self.source_properties = source_properties_
        self.target_properties = target_properties_
        self.distance_threshold = distance_threshold_
        self.n_identity = 0
        self.n_different = 0
        self.CA = 0
        self.CZ = 0
        self.N = 0
        self.NZ = 0
        self.O = 0
        self.OG = 0
        self.OD1 = 0
        self.DU = 0

        if len(self.source_coordinates) > len(self.target_coordinates):
            N = NearestNeighbors(n_neighbors=len(self.source_coordinates),
                      radius = self.distance_threshold,
                      algorithm='ball_tree').fit(self.source_coordinates)
            self.distances, self.indices = N.radius_neighbors(self.target_coordinates)
            self.fitProp = self.target_properties
            self.refProp = self.source_properties

        else:
            N = NearestNeighbors(n_neighbors=len(self.target_coordinates),
                                 radius = self.distance_threshold,
                                 algorithm='ball_tree').fit(self.target_coordinates)
            self.distances, self.indices = N.radius_neighbors(self.source_coordinates)
            self.fitProp = self.source_properties
            self.refProp = self.target_properties

        self.fitsize = len(self.fitProp)
        self.refsize = len(self.refProp)


    def get_similarity_by_rules(self):
        
        for i in range(len(self.fitProp)):

            if self.fitProp[i][1] in [self.refProp[e][1] for e in self.indices[i]]:
                #print(self.fitProp[i][1], self.distances[i])
                """print(self.fitProp[i][1], 
                     [self.refProp[e][1] for e in self.indices[i]].count(self.fitProp[i][1]),

                     [self.refProp[self.indices[i][e]][1] 
                        for e in range(len(self.indices[i]))
                        if self.distances[i][e] <= args.distance].count("CA"),
                     [self.refProp[self.indices[i][e]][1] 
                        for e in range(len(self.indices[i]))
                        if self.distances[i][e] <= args.distance].count("CZ"),
                     [self.refProp[self.indices[i][e]][1] 
                        for e in range(len(self.indices[i]))
                        if self.distances[i][e] <= args.distance].count("N"),
                     [self.refProp[self.indices[i][e]][1] 
                        for e in range(len(self.indices[i]))
                        if self.distances[i][e] <= args.distance].count("NZ"),
                     [self.refProp[self.indices[i][e]][1] 
                        for e in range(len(self.indices[i]))
                        if self.distances[i][e] <= args.distance].count("O"),
                     [self.refProp[self.indices[i][e]][1] 
                        for e in range(len(self.indices[i]))
                        if self.distances[i][e] <= args.distance].count("OD1"),
                     [self.refProp[self.indices[i][e]][1] 
                        for e in range(len(self.indices[i]))
                        if self.distances[i][e] <= args.distance].count("OG"),
                     [self.refProp[self.indices[i][e]][1] 
                        for e in range(len(self.indices[i]))
                        if self.distances[i][e] <= args.distance].count("DU"),
                      self.distances[i]

                     )"""
                self.n_identity += 1
                if self.refProp[i][1] == "CA":
                    self.CA += 1
                if self.refProp[i][1] == "CZ":
                    self.CZ += 1
                if self.refProp[i][1] == "N":
                    self.N += 1
                if self.refProp[i][1] == "NZ":
                    self.NZ += 1
                if self.refProp[i][1] == "O":
                    self.O += 1
                if self.refProp[i][1] == "OD1":
                    self.OD1 += 1
                if self.refProp[i][1] == "OG":
                    self.OG += 1
                if self.refProp[i][1] == "DU":
                    self.DU += 1
                    
        ratio_aligned = round(float(self.n_identity)/self.fitsize, 4)
        
        if self.n_identity != 0:
            ratio_CA_in_aligned = round(float(self.CA)/self.n_identity, 4)
            ratio_CZ_in_aligned = round(float(self.CZ)/self.n_identity, 4)
            ratio_N_in_aligned = round(float(self.N)/self.n_identity, 4)
            ratio_NZ_in_aligned = round(float(self.NZ)/self.n_identity, 4)
            ratio_O_in_aligned = round(float(self.O)/self.n_identity, 4)
            ratio_OD1_in_aligned = round(float(self.OD1)/self.n_identity, 4)
            ratio_OG_in_aligned = round(float(self.OG)/self.n_identity, 4)
            ratio_DU_in_aligned = round(float(self.DU)/self.n_identity, 4)
        else:
            ratio_CA_in_aligned = 0
            ratio_CZ_in_aligned = 0
            ratio_N_in_aligned = 0
            ratio_NZ_in_aligned = 0
            ratio_O_in_aligned = 0
            ratio_OD1_in_aligned = 0
            ratio_OG_in_aligned = 0
            ratio_DU_in_aligned = 0

        return ratio_aligned, \
            ratio_CA_in_aligned, ratio_CZ_in_aligned, \
            ratio_N_in_aligned, ratio_NZ_in_aligned, \
            ratio_O_in_aligned, ratio_OD1_in_aligned, \
            ratio_OG_in_aligned, ratio_DU_in_aligned





class _ph4_soft_(_similarity_metrics_, _distances_):
    """ 1-NN search and 
        strict correspondence of properties """


    def __init__(self, source_coordinates_, target_coordinates_, 
                 source_properties_, target_properties_, distance_threshold_):

        import numpy as np
        from sklearn.neighbors import NearestNeighbors
        
        self.source_coordinates = source_coordinates_
        self.target_coordinates = target_coordinates_
        self.source_properties = source_properties_
        self.target_properties = target_properties_
        self.distance_threshold = distance_threshold_
        self.n_identity = 0
        self.n_different = 0
        self.CA = 0
        self.CZ = 0
        self.N = 0
        self.NZ = 0
        self.O = 0
        self.OG = 0
        self.OD1 = 0
        self.DU = 0

        if len(self.source_coordinates) > len(self.target_coordinates):
            N = NearestNeighbors(n_neighbors=1, 
                                 algorithm="ball_tree").fit(self.source_coordinates)
            self.distances, self.indices = N.kneighbors(np.array(self.target_coordinates))
            self.fitProp = self.target_properties
            self.refProp = self.source_properties
        
        else:
            N = NearestNeighbors(n_neighbors=1, 
                                 algorithm="ball_tree").fit(self.target_coordinates)
            self.distances, self.indices = N.kneighbors(np.array(self.source_coordinates))
            self.fitProp = self.source_properties
            self.refProp = self.target_properties

        self.fitsize = len(self.fitProp)
        self.refsize = len(self.refProp)


    def get_similarity_by_rules(self):
        
        for i in range(len(self.fitProp)):

            if (self.fitProp[i][1] == "CA") and\
                    (self.refProp[self.indices[i][0]][1] == "CA"):
                self.n_identity += 1
                self.CA += 1

            elif (self.fitProp[i][1] == "CZ") and\
                    (self.refProp[self.indices[i][0]][1] == "CZ"):
                self.n_identity += 1
                self.CZ += 1

            elif (self.fitProp[i][1] == "N") and\
                    (self.refProp[self.indices[i][0]][1] == "N"):
                self.n_identity += 1
                self.N += 1
                
            elif (self.fitProp[i][1] == "O") and\
                    (self.refProp[self.indices[i][0]][1] == "O"):
                self.n_identity += 1
                self.O += 1

            elif (self.fitProp[i][1] == "NZ") and\
                    (self.refProp[self.indices[i][0]][1] == "NZ"):
                self.n_identity += 1
                self.NZ += 1

            elif (self.fitProp[i][1] == "OD1") and\
                    (self.refProp[self.indices[i][0]][1] == "OD1"):
                self.n_identity += 1
                self.OD1 += 1

            elif (self.fitProp[i][1] == "OG") and\
                    (self.refProp[self.indices[i][0]][1] == "OG"):
                self.n_identity += 1
                self.OG += 1

            elif (self.fitProp[i][1] == "DU") and\
                    (self.refProp[self.indices[i][0]][1] == "DU"):
                self.n_identity += 1
                self.DU += 1

        ratio_aligned = round(float(self.n_identity)/self.fitsize, 4)
        
        if self.n_identity != 0:
            ratio_CA_in_aligned = round(float(self.CA)/self.n_identity, 4)
            ratio_CZ_in_aligned = round(float(self.CZ)/self.n_identity, 4)
            ratio_N_in_aligned = round(float(self.N)/self.n_identity, 4)
            ratio_NZ_in_aligned = round(float(self.NZ)/self.n_identity, 4)
            ratio_O_in_aligned = round(float(self.O)/self.n_identity, 4)
            ratio_OD1_in_aligned = round(float(self.OD1)/self.n_identity, 4)
            ratio_OG_in_aligned = round(float(self.OG)/self.n_identity, 4)
            ratio_DU_in_aligned = round(float(self.DU)/self.n_identity, 4)
        else:
            ratio_CA_in_aligned = 0
            ratio_CZ_in_aligned = 0
            ratio_N_in_aligned = 0
            ratio_NZ_in_aligned = 0
            ratio_O_in_aligned = 0
            ratio_OD1_in_aligned = 0
            ratio_OG_in_aligned = 0
            ratio_DU_in_aligned = 0

        return ratio_aligned, \
            ratio_CA_in_aligned, ratio_CZ_in_aligned, \
            ratio_N_in_aligned, ratio_NZ_in_aligned, \
            ratio_O_in_aligned, ratio_OD1_in_aligned, \
            ratio_OG_in_aligned, ratio_DU_in_aligned



if __name__ == '__main__':

    import argparse


