

__author__ = "Merveille Eguida"
__copyright__ = "2020"
__license__ = "MIT"
__version__ = "0.1.2"
__maintainer__ = "Merveille Eguida"
__email__ = "keguida@unistra.fr"
__status__ = "Development"



import os
import argparse
import numpy as np
import math
from sklearn.neighbors import NearestNeighbors





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
                if self.fitProp[i][1] == "CA":
                    self.CA += 1
                if self.fitProp[i][1] == "CZ":
                    self.CZ += 1
                if self.fitProp[i][1] == "N":
                    self.N += 1
                if self.fitProp[i][1] == "NZ":
                    self.NZ += 1
                if self.fitProp[i][1] == "O":
                    self.O += 1
                if self.fitProp[i][1] == "OD1":
                    self.OD1 += 1
                if self.fitProp[i][1] == "OG":
                    self.OG += 1
                if self.fitProp[i][1] == "DU":
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




class _piecewise_linear_:

    def piecewise_linear_func(self, distance_, distance_threshold_):
        if 0 <= distance_ < float(distance_threshold_)/2:
            pl = 1
        elif float(distance_threshold_)/2 <= distance_ < distance_threshold_:
            pl = 1 - float(distance_)/distance_threshold_
        elif distance_ >= distance_threshold_:
            pl = 0

        return pl




class _ph4_strict_pl_(_ph4_strict_, _piecewise_linear_):
    # PL: piecewise linear

    def get_similarity_by_rules(self):
        
        for i in range(len(self.fitProp)):

            if (self.fitProp[i][1] == "CA") and\
                    (self.refProp[self.indices[i][0]][1] == "CA"):
                self.n_identity += self.piecewise_linear_func(
                        self.distances[i][0], self.distance_threshold)
                self.CA += 1

            elif (self.fitProp[i][1] == "CZ") and\
                    (self.refProp[self.indices[i][0]][1] == "CZ"):
                self.n_identity += self.piecewise_linear_func(
                        self.distances[i][0], self.distance_threshold)
                self.CZ += 1

            elif (self.fitProp[i][1] == "N") and\
                    (self.refProp[self.indices[i][0]][1] == "N"):
                self.n_identity += self.piecewise_linear_func(
                        self.distances[i][0], self.distance_threshold)
                self.N += 1
                
            elif (self.fitProp[i][1] == "O") and\
                    (self.refProp[self.indices[i][0]][1] == "O"):
                self.n_identity += self.piecewise_linear_func(
                        self.distances[i][0], self.distance_threshold)
                self.O += 1

            elif (self.fitProp[i][1] == "NZ") and\
                    (self.refProp[self.indices[i][0]][1] == "NZ"):
                self.n_identity += self.piecewise_linear_func(
                        self.distances[i][0], self.distance_threshold)
                self.NZ += 1

            elif (self.fitProp[i][1] == "OD1") and\
                    (self.refProp[self.indices[i][0]][1] == "OD1"):
                self.n_identity += self.piecewise_linear_func(
                        self.distances[i][0], self.distance_threshold)
                self.OD1 += 1

            elif (self.fitProp[i][1] == "OG") and\
                    (self.refProp[self.indices[i][0]][1] == "OG"):
                self.n_identity += self.piecewise_linear_func(
                        self.distances[i][0], self.distance_threshold)
                self.OG += 1

            elif (self.fitProp[i][1] == "DU") and\
                    (self.refProp[self.indices[i][0]][1] == "DU"):
                self.n_identity += self.piecewise_linear_func(
                        self.distances[i][0], self.distance_threshold)
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
    



class _ph4_rules_pl_(_ph4_rules_, _piecewise_linear_):

    def get_similarity_by_rules(self):
        
        for i in range(len(self.fitProp)):
            
            if (self.fitProp[i][1] == "CA") and\
                    (self.refProp[self.indices[i][0]][1] == "CA" or "CZ"):
                self.n_identity += self.piecewise_linear_func(
                        self.distances[i][0], self.distance_threshold)
                self.CA += 1

            elif (self.fitProp[i][1] == "CZ") and\
                    (self.refProp[self.indices[i][0]][1] == "CA" or "CZ"):
                self.n_identity += self.piecewise_linear_func(
                        self.distances[i][0], self.distance_threshold)
                self.CZ += 1
 
            elif (self.fitProp[i][1] == "N") and\
                    (self.refProp[self.indices[i][0]][1] == "N" or "NZ" or "OG"):
                self.n_identity += self.piecewise_linear_func(
                        self.distances[i][0], self.distance_threshold)
                self.N += 1
                
            elif (self.fitProp[i][1] == "O") and\
                    (self.refProp[self.indices[i][0]][1] == "O" or "OD1" or "OG"):
                self.n_identity += self.piecewise_linear_func(
                        self.distances[i][0], self.distance_threshold)
                self.O += 1

            elif (self.fitProp[i][1] == "NZ") and\
                    (self.refProp[self.indices[i][0]][1] == "NZ" or "N" or "OG"):
                self.n_identity += self.piecewise_linear_func(
                        self.distances[i][0], self.distance_threshold)
                self.NZ += 1

            elif (self.fitProp[i][1] == "OD1") and\
                    (self.refProp[self.indices[i][0]][1] == "OD1" or "O" or "OG"):
                self.n_identity += self.piecewise_linear_func(
                        self.distances[i][0], self.distance_threshold)
                self.OD1 += 1

            elif (self.fitProp[i][1] == "OG") and\
                    (self.refProp[self.indices[i][0]][1] == "OG" or "N" or "O" or "NZ" or "OD1"):
                self.n_identity += self.piecewise_linear_func(
                        self.distances[i][0], self.distance_threshold)
                self.OG += 1

            elif (self.fitProp[i][1] == "DU") and\
                    (self.refProp[self.indices[i][0]][1] == "DU"):
                self.n_identity += self.piecewise_linear_func(
                        self.distances[i][0], self.distance_threshold)
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




class _ph4_ext_pl_(_ph4_ext_, _piecewise_linear_):

    def get_similarity_by_rules(self):
        
        for i in range(len(self.fitProp)):

            if self.fitProp[i][1] in [self.refProp[e][1] for e in self.indices[i]]:
                self.n_identity += self.piecewise_linear_func(
                        self.distances[i][0], self.distance_threshold)
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




class _fingerprint_distances_rules_:

    def __search_1NN(self, ref_vector_, fit_vector_):
        N = NearestNeighbors(n_neighbors=1, 
                                algorithm="ball_tree").fit(ref_vector_)
        distances, indices = N.kneighbors(np.array(fit_vector_))
        np.reshape(distances, len(fit_vector_))
        return distances


    def distance_rules(self, method_, ref_vector_, fit_vector_):
        """ # all: average of the 1-NN distances
            # q3: average of the 1-NN distances inferior to the 3rd quantile Q3
            # q2 average of the 1-NN distances inferior to median Q2
        """
        if method_ == 'all':
            distances = self.__search_1NN(ref_vector_, fit_vector_)
            dist = np.mean(distances)


        elif method_ == 'q3':
            distances = self.__search_1NN(ref_vector_, fit_vector_)
            q3_distances = [d for d in distances 
                            if d >= np.quantile(distances, 0.75)]
            dist = np.mean(q3_distances)


        elif method_ == 'q2':
            distances = self.__search_1NN(ref_vector_, fit_vector_)
            q2_distances = [d for d in distances 
                            if d >= np.quantile(distances, 0.50)]
            dist = np.mean(q2_distances)

        return dist



class _color_fingerprint_distances_(_fingerprint_distances_rules_):
    """ """

    def __init__(self, source_coordinates_, target_coordinates_, 
                 source_properties_, target_properties_, distance_threshold_):
        self.source_coordinates = source_coordinates_
        self.target_coordinates = target_coordinates_
        self.source_properties = source_properties_
        self.target_properties = target_properties_
        self.distance_threshold = distance_threshold_

        if len(self.source_coordinates) > len(self.target_coordinates):
            self.ref_coordinates = self.source_coordinates
            self.ref_properties = self.source_properties
            self.fit_coordinates = self.target_coordinates
            self.fit_properties = self.target_properties
        else:
            self.ref_coordinates = self.target_coordinates
            self.ref_properties = self.target_properties
            self.fit_coordinates = self.source_coordinates
            self.fit_properties = self.source_properties



    def __fingerprint(self, coordinates_, properties_):
        
        N = NearestNeighbors(n_neighbors=len(coordinates_), # max
                      radius = self.distance_threshold,
                      algorithm='ball_tree').fit(coordinates_)

        distances, indices = N.radius_neighbors(coordinates_)
        fingerprint = []
        for i in range(len(coordinates_)):
            fp = [0 for i in range(8)]
            for neig_indx in indices[i]:
                if properties_[neig_indx][1] == "CA":
                    fp[0] += 1
                elif properties_[neig_indx][1] == "CZ":
                    fp[1] += 1
                elif properties_[neig_indx][1] == "O":
                    fp[2] += 1
                elif properties_[neig_indx][1] == "OD1":
                    fp[3] += 1
                elif properties_[neig_indx][1] == "OG":
                    fp[4] += 1
                elif properties_[neig_indx][1] == "N":
                    fp[5] += 1
                elif properties_[neig_indx][1] == "NZ":
                    fp[6] += 1
                elif properties_[neig_indx][1] == "DU":
                    fp[7] += 1

            fp = np.array(fp) / len(indices[i]) * 100
            fingerprint.append(fp)
        if np.sum(fingerprint) == 0:
            print(fingerprint)
        return np.array(fingerprint)




    def get_distance_by_rules(self, method_='all'):
        
        self.ref_fingerprint = self.__fingerprint(self.ref_coordinates, 
                                                  self.ref_properties)
        self.fit_fingerprint = self.__fingerprint(self.fit_coordinates, 
                                                  self.fit_properties)

        dist = self.distance_rules(method_=method_,
                                    ref_vector_=self.ref_fingerprint,
                                    fit_vector_=self.fit_fingerprint)
        
        self.distance = dist


def read_cav_mol2(mol2_):

    with open(mol2_, "r") as mol2file:
        atoms = mol2file.read().split("@<TRIPOS>ATOM\n")[-1].split('@<TRIPOS>')[0]
        atoms = atoms.split('\n')
        del atoms[-1]
        
        crd = []
        properties = []
        for atm_line in atoms:
            if atm_line == '':
                continue
            cols = atm_line.split()
            idx = int(cols[0])
            atom_type = str(cols[1])
            x = float(cols[2])
            y = float(cols[3])
            z = float(cols[4])
            crd.append([x,y,z])
            properties.append([idx, atom_type])

    return np.array(crd), properties


if __name__ == '__main__':


    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--source', type=str, required=True, 
                         help="source / input cavity mol2 (1)")
    parser.add_argument('-t', '--target', type=str, required=True, 
                         help="target / input cavity mol2 (2)")
    parser.add_argument('-d', '--distance', type=float, required=False, 
                         default=1.5, help="Distance threshold")
    parser.add_argument('-o', '--ofile', type=str, required=False,
                         default='procare_rescoring.tsv', help='rescoring output file')
    args = parser.parse_args()


    source_coordinates, source_properties = read_cav_mol2(args.source)
    target_coordinates, target_properties = read_cav_mol2(args.target)



    ph4_strict = _ph4_strict_(source_coordinates, target_coordinates, 
            source_properties, target_properties, args.distance)
    ratio_aligned, ratio_CA_in_aligned, \
        ratio_CZ_in_aligned, ratio_N_in_aligned, \
        ratio_NZ_in_aligned, ratio_O_in_aligned, \
        ratio_OD1_in_aligned, ratio_OG_in_aligned, \
        ratio_DU_in_aligned = ph4_strict.get_similarity_by_rules()


    ph4_rules = _ph4_rules_(source_coordinates, target_coordinates, 
            source_properties, target_properties, args.distance)
    ratio_aligned, ratio_CA_in_aligned, \
        ratio_CZ_in_aligned, ratio_N_in_aligned, \
        ratio_NZ_in_aligned, ratio_O_in_aligned, \
        ratio_OD1_in_aligned, ratio_OG_in_aligned, \
        ratio_DU_in_aligned = ph4_rules.get_similarity_by_rules()


    ph4_ext = _ph4_ext_(source_coordinates, target_coordinates, 
            source_properties, target_properties, args.distance)
    ratio_aligned, ratio_CA_in_aligned, \
        ratio_CZ_in_aligned, ratio_N_in_aligned, \
        ratio_NZ_in_aligned, ratio_O_in_aligned, \
        ratio_OD1_in_aligned, ratio_OG_in_aligned, \
        ratio_DU_in_aligned = ph4_ext.get_similarity_by_rules()


    ph4_soft = _ph4_soft_(source_coordinates, target_coordinates, 
            source_properties, target_properties, args.distance)
    ratio_aligned, ratio_CA_in_aligned, \
        ratio_CZ_in_aligned, ratio_N_in_aligned, \
        ratio_NZ_in_aligned, ratio_O_in_aligned, \
        ratio_OD1_in_aligned, ratio_OG_in_aligned, \
        ratio_DU_in_aligned = ph4_soft.get_similarity_by_rules()


    source_filename = os.path.basename(args.source)
    target_filename = os.path.basename(args.target)

    if not os.path.isfile(args.ofile):
        with open(args.ofile, 'w') as of:
            of.write('source\ttarget\t'
                'ph4_strict_tv\tph4_ext_tv\tph4_soft_tv\tph4_rules_tv\t'
                'ph4_strict_tc\tph4_ext_tc\tph4_soft_tc\tph4_rules_tc\n')

    with open(args.ofile, 'a') as of:
        of.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
            source_filename,
            target_filename,
            ph4_strict.tversky_similarity(),
            ph4_ext.tversky_similarity(),
            ph4_soft.tversky_similarity(),
            ph4_rules.tversky_similarity(),
            ph4_strict.tanimoto_similarity(),
            ph4_ext.tanimoto_similarity(),
            ph4_soft.tanimoto_similarity(),
            ph4_rules.tanimoto_similarity()
            ))

    """
    print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                source_filename,
                target_filename,
                ph4_strict.tversky_similarity(),
                ph4_ext.tversky_similarity(),
                ph4_soft.tversky_similarity(),
                ph4_rules.tversky_similarity(),
                ph4_strict.tanimoto_similarity(),
                ph4_ext.tanimoto_similarity(),
                ph4_soft.tanimoto_similarity(),
                ph4_rules.tanimoto_similarity()
                ))
    """