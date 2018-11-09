#ifndef LAB4INTERSECTION
#define LAB4INTERSECTION

#include <iostream>
#include <Eigen\Dense> //for Eigen library
#include <Eigen\Core>
#include <fstream> //for read/write files
#include <stdlib.h>// for system("pause") and exit_failure
#include <vector>
#include <math.h>  
#include <algorithm>
#include <stdexcept>
#include <iomanip> 
#include <string>


using namespace std;
using namespace Eigen;

#define MaxMatSize 10000000
#define pi 3.14159265358979323846 
struct CameraParam {
	double PS;
	double f_l;
	double xpp ;
	double ypp;
	double K1 ;
	double K2 ;
	double K3 ;
	double P1 ;
	double P2 ;
	double S1 ;
	double S2 ;
	double Cn ;
	double Rn ;

};

typedef Matrix<double, 2, 1> Matrix2b1;
typedef Matrix<double, Dynamic, 3> Matrixdb3;
typedef Matrix<double, 3, 3> Matrix3b3;
typedef Matrix<double, 3, 4> Matrix3b4;
typedef Matrix<double, Dynamic, 2> Matrixdby2;
typedef Matrix<int, Dynamic, 2> Matrixdby2i;

// Read and write function are modified versions of DMP_BBO library (TAKEN FROM DR.SHAHBAZI)
void Read_Mat(char *FileName, MatrixXd& m);//reads a formatted file to a matrix
void Write_Mat(char *FileName, MatrixXd & m, int decimal_precision); //writes a matrix to a formatted file

void SplitObs_and_RANSAC(CameraParam camera_params, MatrixXd tie_pts, double ransac_conf, double outlier_percentage, double min_iterations, double dThreshold);

void Rotation_g2i(double Omega, double Phi, double Kappa, Matrix3b3 & Rot_g2i);

void removeRow(MatrixXd& matrix, unsigned int rowToRemove);//remove a row from a matrix
void removeColumn(MatrixXd& matrix, unsigned int colToRemove);//remove a column from a matrix

void  Normalization_Condition(MatrixXd &xy_i1, MatrixXd &xy_i2, MatrixXd& H1, MatrixXd& H2);

void Perform_LinOri(CameraParam& camera_params, MatrixXd &xy_i1, MatrixXd &xy_i2, MatrixXd& Fun, MatrixXd& Emat);

bool Decompose_Essential(CameraParam& camera_params, MatrixXd Emat, MatrixXd xy_i1, MatrixXd xy_i2,	double& Bx, double& By, double& Bz, double& Omega, double& Phi, double& Kappa);

bool Perform_NonlinOri(CameraParam& camera_params, MatrixXd &xy_i1, MatrixXd &xy_i2, MatrixXd& T21, double& Bx, double& By, double& Bz, double& Omega, double& Phi, double& Kappa);

void Vanilla_RANSAC(CameraParam& camera_params, MatrixXd &xy_i1, MatrixXd &xy_i2, double confidence_level, double outlier_percentage, int min_iterations, double dThreshold, MatrixXd& Inliers_index);

void FindInliers(MatrixXd inliers, MatrixXd all_ties_1, MatrixXd all_ties_2, MatrixXd &inlier_ties1, MatrixXd &inlier_ties2);

#endif
