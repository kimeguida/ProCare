// ########################## OPEN3D ORIGINAL WORK ############################
// ----------------------------------------------------------------------------
// -                        Open3D: www.open3d.org                            -
// ----------------------------------------------------------------------------
// The MIT License (MIT)
//
// Copyright (c) 2018 www.open3d.org
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.
// ----------------------------------------------------------------------------
//
//
//
// ########################## PROCARE MODIFIED WORK ############################
// Mofifications: block of codes specified by /*kimeguida*/
// -----------------------------------------------------------------------------
// <                                  ProCare                                  >
// -----------------------------------------------------------------------------
// The MIT License (MIT)
//
// Copyright (c) 2020 Merveille Eguida
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.


#include "Feature.h"

#include <Eigen/Dense>
#include <Core/Utility/Console.h>
#include <Core/Geometry/PointCloud.h>
#include <Core/Geometry/KDTreeFlann.h>

namespace open3d {

namespace {

Eigen::Vector4d ComputePairFeatures(const Eigen::Vector3d &p1,
                                    const Eigen::Vector3d &n1,
                                    const Eigen::Vector3d &p2,
                                    const Eigen::Vector3d &n2) {
    Eigen::Vector4d result;
    Eigen::Vector3d dp2p1 = p2 - p1;
    result(3) = dp2p1.norm();
    if (result(3) == 0.0) {
        return Eigen::Vector4d::Zero();
    }
    auto n1_copy = n1;
    auto n2_copy = n2;
    double angle1 = n1_copy.dot(dp2p1) / result(3);
    double angle2 = n2_copy.dot(dp2p1) / result(3);
    if (acos(fabs(angle1)) > acos(fabs(angle2))) {
        n1_copy = n2;
        n2_copy = n1;
        dp2p1 *= -1.0;
        result(2) = -angle2;
    } else {
        result(2) = angle1;
    }
    auto v = dp2p1.cross(n1_copy);
    double v_norm = v.norm();
    if (v_norm == 0.0) {
        return Eigen::Vector4d::Zero();
    }
    v /= v_norm;
    auto w = n1_copy.cross(v);
    result(1) = v.dot(n2_copy);
    result(0) = atan2(w.dot(n2_copy), n1_copy.dot(n2_copy));
    return result;
}

std::shared_ptr<Feature> ComputeSPFHFeature(
        const PointCloud &input,
        const KDTreeFlann &kdtree,
        const KDTreeSearchParam &search_param) {
    auto feature = std::make_shared<Feature>();
    feature->Resize(33, (int)input.points_.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < (int)input.points_.size(); i++) {
        const auto &point = input.points_[i];
        const auto &normal = input.normals_[i];
        std::vector<int> indices;
        std::vector<double> distance2;
        if (kdtree.Search(point, search_param, indices, distance2) > 1) {
            // only compute SPFH feature when a point has neighbors
            double hist_incr = 100.0 / (double)(indices.size() - 1);
            for (size_t k = 1; k < indices.size(); k++) {
                // skip the point itself, compute histogram
                auto pf = ComputePairFeatures(point, normal,
                                              input.points_[indices[k]],
                                              input.normals_[indices[k]]);
                int h_index = (int)(floor(11 * (pf(0) + M_PI) / (2.0 * M_PI)));
                if (h_index < 0) h_index = 0;
                if (h_index >= 11) h_index = 10;
                feature->data_(h_index, i) += hist_incr;
                h_index = (int)(floor(11 * (pf(1) + 1.0) * 0.5));
                if (h_index < 0) h_index = 0;
                if (h_index >= 11) h_index = 10;
                feature->data_(h_index + 11, i) += hist_incr;
                h_index = (int)(floor(11 * (pf(2) + 1.0) * 0.5));
                if (h_index < 0) h_index = 0;
                if (h_index >= 11) h_index = 10;
                feature->data_(h_index + 22, i) += hist_incr;
            }
        }
    }
    return feature;
}

/*kimeguida*/
std::shared_ptr<Feature> ComputeCSPFHFeature(
        const PointCloud &input,
        const KDTreeFlann &kdtree,
        const KDTreeSearchParam &search_param) {
    auto feature = std::make_shared<Feature>();
    feature->Resize(41, (int)input.points_.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < (int)input.points_.size(); i++) {
        const auto &point = input.points_[i];
        const auto &normal = input.normals_[i];
        std::vector<int> indices;
        std::vector<double> distance2;
        if (kdtree.Search(point, search_param, indices, distance2) > 1) {
            // only compute CSPFH feature when a point has neighbors
            double hist_incr = 100.0 / (double)(indices.size() - 1);
            for (size_t k = 1; k < indices.size(); k++) {
                // skip the point itself, compute histogram
                auto pf = ComputePairFeatures(point, normal,
                                              input.points_[indices[k]],
                                              input.normals_[indices[k]]);
                int h_index = (int)(floor(11 * (pf(0) + M_PI) / (2.0 * M_PI)));
                if (h_index < 0) h_index = 0;
                if (h_index >= 11) h_index = 10;
                feature->data_(h_index, i) += hist_incr;
                h_index = (int)(floor(11 * (pf(1) + 1.0) * 0.5));
                if (h_index < 0) h_index = 0;
                if (h_index >= 11) h_index = 10;
                feature->data_(h_index + 11, i) += hist_incr;
                h_index = (int)(floor(11 * (pf(2) + 1.0) * 0.5));
                if (h_index < 0) h_index = 0;
                if (h_index >= 11) h_index = 10;
                feature->data_(h_index + 22, i) += hist_incr;
            }

            double c_hist_incr = 100.0 / (double)(indices.size());
            const char *CA_c = "16741671";
            const char *CZ_c = "4646984";
            const char *O_c = "15219528";
            const char *OD1_c = "0";
            const char *OG_c = "8204959";
            const char *N_c = "30894";
            const char *NZ_c = "15231913";
            const char *DU_c = "7566712";

            Eigen::Vector3d CA_rgb = ASCIIPCDColorToRGB(CA_c, 'F', 4);
            Eigen::Vector3d CZ_rgb = ASCIIPCDColorToRGB(CZ_c, 'F', 4);
            Eigen::Vector3d O_rgb = ASCIIPCDColorToRGB(O_c, 'F', 4);
            Eigen::Vector3d OD1_rgb = ASCIIPCDColorToRGB(OD1_c, 'F', 4);
            Eigen::Vector3d OG_rgb = ASCIIPCDColorToRGB(OG_c, 'F', 4);
            Eigen::Vector3d N_rgb = ASCIIPCDColorToRGB(N_c, 'F', 4);
            Eigen::Vector3d NZ_rgb = ASCIIPCDColorToRGB(NZ_c, 'F', 4);
            Eigen::Vector3d DU_rgb = ASCIIPCDColorToRGB(DU_c, 'F', 4);

            for (size_t k = 0; k < indices.size(); k++) {
                // include the point itself, compute histogram
                if (input.colors_[indices[k]] == CA_rgb) {
                    int c_h_index = 33;
                    feature->data_(c_h_index, i) += c_hist_incr;
                }
                if (input.colors_[indices[k]] == CZ_rgb) {
                    int c_h_index = 34;
                    feature->data_(c_h_index, i) += c_hist_incr;
                }
                if (input.colors_[indices[k]] == O_rgb) {
                    int c_h_index = 35;
                    feature->data_(c_h_index, i) += c_hist_incr;
                }
                if (input.colors_[indices[k]] == OD1_rgb) {
                    int c_h_index = 36;
                    feature->data_(c_h_index, i) += c_hist_incr;
                }
                if (input.colors_[indices[k]] == OG_rgb) {
                    int c_h_index = 37;
                    feature->data_(c_h_index, i) += c_hist_incr;
                }
                if (input.colors_[indices[k]] == N_rgb) {
                    int c_h_index = 38;
                    feature->data_(c_h_index, i) += c_hist_incr;
                }
                if (input.colors_[indices[k]] == NZ_rgb) {
                    int c_h_index = 39;
                    feature->data_(c_h_index, i) += c_hist_incr;
                }
                if (input.colors_[indices[k]] == DU_rgb) {
                    int c_h_index = 40;
                    feature->data_(c_h_index, i) += c_hist_incr;
                }
            }
        }
    }
    return feature;
}
/*kimeguida*/

}  // unnamed namespace

std::shared_ptr<Feature> ComputeFPFHFeature(
        const PointCloud &input,
        const KDTreeSearchParam &search_param /* = KDTreeSearchParamKNN()*/) {
    auto feature = std::make_shared<Feature>();
    feature->Resize(33, (int)input.points_.size());
    if (input.HasNormals() == false) {
        PrintDebug(
                "[ComputeFPFHFeature] Failed because input point cloud has no "
                "normal.\n");
        return feature;
    }
    KDTreeFlann kdtree(input);
    auto spfh = ComputeSPFHFeature(input, kdtree, search_param);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < (int)input.points_.size(); i++) {
        const auto &point = input.points_[i];
        std::vector<int> indices;
        std::vector<double> distance2;
        if (kdtree.Search(point, search_param, indices, distance2) > 1) {
            double sum[3] = {0.0, 0.0, 0.0};
            for (size_t k = 1; k < indices.size(); k++) {
                // skip the point itself
                double dist = distance2[k];
                if (dist == 0.0) continue;
                for (int j = 0; j < 33; j++) {
                    double val = spfh->data_(j, indices[k]) / dist;
                    sum[j / 11] += val;
                    feature->data_(j, i) += val;
                }
            }
            for (int j = 0; j < 3; j++)
                if (sum[j] != 0.0) sum[j] = 100.0 / sum[j];
            for (int j = 0; j < 33; j++) {
                feature->data_(j, i) *= sum[j / 11];
                // The commented line is the fpfh function in the paper.
                // But according to PCL implementation, it is skipped.
                // Our initial test shows that the full fpfh function in the
                // paper seems to be better than PCL implementation. Further
                // test required.
                feature->data_(j, i) += spfh->data_(j, i);
            }
        }
    }
    return feature;
}

/*kimeguida*/
Eigen::Vector3d ASCIIPCDColorToRGB(const char *color_ptr,
                                   const char type,
                                   const int size) {
    if ((size == 4) && (type == 'F')) {
        std::uint8_t c_data[4] = {0, 0, 0, 0};
        char *end;
        std::float_t c_value = std::strtof(color_ptr, &end);
        memcpy(c_data, &c_value, 4);
        return Eigen::Vector3d((double)c_data[2] / 255.0,
                               (double)c_data[1] / 255.0,
                               (double)c_data[0] / 255.0);
    } else {
        return Eigen::Vector3d::Zero();
    }
}

std::shared_ptr<Feature> ComputeCFPFHFeature(
        const PointCloud &input,
        const KDTreeSearchParam
                &search_param /* = geometry::KDTreeSearchParamKNN()*/) {
    auto feature = std::make_shared<Feature>();
    feature->Resize(41, (int)input.points_.size());
    if (input.HasNormals() == false) {
        PrintDebug(
                "[ComputeCFPFHFeature] Failed because input point cloud has no "
                "normal.\n");
        return feature;
    }
    if (input.HasColors() == false) {
        PrintDebug(
                "[ComputeCFPFHFeature] Failed because input point cloud has no "
                "color.\n");
        return feature;
    }
    KDTreeFlann kdtree(input);
    auto cspfh = ComputeCSPFHFeature(input, kdtree, search_param);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < (int)input.points_.size(); i++) {
        const auto &point = input.points_[i];
        std::vector<int> indices;
        std::vector<double> distance2;
        if (kdtree.Search(point, search_param, indices, distance2) > 1) {
            double sum[3] = {0.0, 0.0, 0.0};
            double c_sum = 0;
            for (size_t k = 1; k < indices.size(); k++) {
                // skip the point itself
                double dist = distance2[k];
                if (dist == 0.0) continue;
                for (int j = 0; j < 33; j++) {
                    double val = cspfh->data_(j, indices[k]) / dist;
                    sum[j / 11] += val;
                    feature->data_(j, i) += val;
                }
                for (int j = 33; j < 41; j++) {
                    double val = cspfh->data_(j, indices[k]) / dist;
                    c_sum += val;
                    feature->data_(j, i) += val;
                }
            }
            for (int j = 0; j < 3; j++)
                if (sum[j] != 0.0) sum[j] = 100.0 / sum[j];
            for (int j = 0; j < 33; j++) {
                feature->data_(j, i) *= sum[j / 11];
                // The commented line is the fpfh function in the paper.
                // But according to PCL implementation, it is skipped.
                // Our initial test shows that the full fpfh function in the
                // paper seems to be better than PCL implementation. Further
                // test required.
                feature->data_(j, i) += cspfh->data_(j, i);
            }

            if (c_sum != 0.0) c_sum = 100.0 / c_sum;

            for (int j = 33; j < 41; j++) {
                feature->data_(j, i) *= c_sum;
                feature->data_(j, i) += cspfh->data_(j, i);
            }
        }
    }
    return feature;
}
/*kimeguida*/
}  // namespace open3d
