#ifndef CSVTOEH
#define CSVTOEH

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "Types.h"
#include "Constraint.h"

inline void IsOpen_HandleError(bool open, const char* file, int line)
{
    if (!open) {
        std::cout << "File Opening Problem in " << file << " at line " << line
                  << std::endl;
        // printf( "File Opening Problem in %s at line %d\n", file, line );
        exit(EXIT_FAILURE);
    }
}
#define ISOPEN_HANDLE_ERROR(err) (IsOpen_HandleError(err, __FILE__, __LINE__))

static plb::npfem::Matrix3X setPointsFromCSV(std::string filename){

	std::vector<std::vector<double>> p_csv;
    int n_points = 0;
    std::ifstream points_file(filename.c_str());
	//std::cout << filename.c_str() << std::endl;
    ISOPEN_HANDLE_ERROR(points_file.is_open());
    // read every point line-by-line
    std::string line;
    while (std::getline(points_file, line)) {
        std::vector<double> point;
        std::istringstream ss(line);
        std::string token;
        while (std::getline(ss, token, ',')) {
            point.push_back(std::stof(token));
        }
        p_csv.push_back(point);
        n_points++;
    }
    points_file.close();

    plb::npfem::Matrix3X q(3, n_points);
    for (int i = 0; i < q.cols(); ++i) q.col(i) = plb::npfem::Vector3(p_csv[i][0], p_csv[i][1], p_csv[i][2]);

    return q;
}

////////////////////////////////////////////////////////////////////////////////
static std::vector<std::shared_ptr<plb::npfem::Constraint>> setConstraintsFromCSV(std::string filename, const plb::npfem::Matrix3X &points)
{
	std::vector<std::shared_ptr<plb::npfem::Constraint>> constraints;

    std::ifstream constraints_file(filename.c_str());
    ISOPEN_HANDLE_ERROR(constraints_file.is_open()); 
    // read every constraint line-by-line
    std::string line;
    while (std::getline(constraints_file, line)) {
        // read details of the constraint
        std::string constraintType;
        std::vector<int> idI;
        plb::npfem::Scalar weight;
        std::vector<plb::npfem::Scalar> scalars;

        std::istringstream ss(line);
        std::string token;
        int field = 0;
        while (std::getline(ss, token, ',')) {
            if (field == 0) {
                constraintType = token;
            } else if (field == 1) {
                std::istringstream ss_tmp(token);
                std::string inds;
                while (std::getline(ss_tmp, inds, ' ')) {
                    idI.push_back(std::stoi(inds));
                }
            } else if (field == 2) {
                weight = std::stof(token);
            } else {
                if (token.length() != 0) {
                    std::istringstream ss_tmp(token);
                    std::string scalar;
                    while (std::getline(ss_tmp, scalar, ' ')) {
                        scalars.push_back(std::stof(scalar));
                    }
                }
            }
            ++field;
        }
        if (constraintType.compare("VolumeMaterial") == 0) {
            auto c = std::make_shared<plb::npfem::VolumeMaterialConstraint>( idI, weight, points);
            c->setRangeMin(scalars[0]);
            c->setRangeMax(scalars[1]);
            c->setMiu(scalars[2]);
            c->setLambda(scalars[3]);
            c->setKappa(scalars[4]);
			constraints.push_back(c);
        } else if (constraintType.compare("VolumeDamping") == 0) {
            auto c = std::make_shared<plb::npfem::VolumeDampingConstraint>( idI, weight, points);
            c->setRangeMin(scalars[0]);
            c->setRangeMax(scalars[1]);
            c->setMiu(scalars[2]);
            c->setLambda(scalars[3]);
            c->setKappa(scalars[4]);
			constraints.push_back(c);
        } else if (constraintType.compare("SurfaceMaterial") == 0) {
            auto c = std::make_shared<plb::npfem::SurfaceMaterialConstraint>( idI, weight, points);
            c->setRangeMin(scalars[0]);
            c->setRangeMax(scalars[1]);
            c->setMiu(scalars[2]);
            c->setLambda(scalars[3]);
            c->setKappa(scalars[4]);
			constraints.push_back(c);
        } else if (constraintType.compare("Volume") == 0) {
            auto c = std::make_shared<plb::npfem::VolumeConstraint>( idI, weight, points);
            c->setRangeMin(scalars[0]);
            c->setRangeMax(scalars[1]);
			constraints.push_back(c);
        } else if (constraintType.compare("Area") == 0) {
            auto c = std::make_shared<plb::npfem::AreaConstraint>( idI, weight, points);
            c->setRangeMin(scalars[0]);
            c->setRangeMax(scalars[1]);
			constraints.push_back(c);
        } else if (constraintType.compare("Bending") == 0) {
            auto c = std::make_shared<plb::npfem::BendingConstraint>(idI, weight, points);
            c->setRangeMin(scalars[0]);
            c->setRangeMax(scalars[1]);
			constraints.push_back(c);
        } else if (constraintType.compare("EdgeStrainLimiting") == 0) {
            auto c = std::make_shared<plb::npfem::EdgeStrainLimitingConstraint>( idI, weight, points);
            c->setRangeMin(scalars[0]);
            c->setRangeMax(scalars[1]);
			constraints.push_back(c);
        } else if (constraintType.compare("TriangleStrainLimiting") == 0) {
            auto c = std::make_shared<plb::npfem::TriangleStrainLimitingConstraint>( idI, weight, points);
            c->setRangeMin(scalars[0]);
            c->setRangeMax(scalars[1]);
			constraints.push_back(c);
        } else {
            std::cout << "Some Constraints are not valid: Check the "
                         "configuration file!"
                      << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    constraints_file.close();
	return constraints;
}

static std::vector<std::array<int, 3>> setConnectivityListFromCSV(std::string filename){

	std::vector<std::array<int, 3>> connectivity_csv;
    int n_triangles = 0;
    std::ifstream connectivity_file(filename.c_str());
    ISOPEN_HANDLE_ERROR(connectivity_file.is_open());
    std::string line;
    while (std::getline(connectivity_file, line)) {
        std::array<int,3> triangle;
        std::istringstream ss(line);
        std::string token;
		int i = 0;
        while (std::getline(ss, token, ',')) {
            triangle[i++] = std::stoi(token);
        }
        connectivity_csv.push_back(triangle);
        n_triangles++;
    }
    connectivity_file.close();

    return connectivity_csv;
}

std::vector<bool> setOnSurfaceParticle(std::string filename){
    std::vector<bool> onSurfaceParticle;

    std::ifstream onSurfaceParticle_file(filename.c_str());
    ISOPEN_HANDLE_ERROR(onSurfaceParticle_file.is_open());

    std::string line;
    while (std::getline(onSurfaceParticle_file, line)) {
        onSurfaceParticle.push_back(std::stoi(line));
    }

    onSurfaceParticle_file.close();

   return onSurfaceParticle;
}

#endif
