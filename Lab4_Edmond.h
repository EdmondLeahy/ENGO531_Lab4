#ifndef LAB4_EDMOND
#define LAB4_EDMOND

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
	double PS; //pixel size
	double f_l; //focal length
	double xpp ; //principal point in x
	double ypp; //principal point in y
	double K1 ; //Distortion parameter
	double K2 ;//Distortion parameter
	double K3 ;//Distortion parameter
	double P1 ;//Distortion parameter
	double P2 ;//Distortion parameter
	double S1 ;//Distortion parameter
	double S2 ;//Distortion parameter
	double Cn ;
	double Rn ;
	double sigma_obs; //observation standard deviation

};
struct RelativeOrientation {
	int from_img;
	int to_img;
	double bx;
	double by;
	double bz;
	double omega;
	double phi;
	double kappa;
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

void SplitObs(CameraParam camera_params, MatrixXd tie_pts, vector<MatrixXd>& all_split_obs, MatrixXd& img_indeces);

void Ransac_All_obs(CameraParam camera_params, vector<MatrixXd>& all_split_obs, double ransac_conf, double outlier_percentage, double min_iterations, double dThreshold);

void Rotation_g2i(double Omega, double Phi, double Kappa, Matrix3b3 & Rot_g2i);

void removeRow(MatrixXd& matrix, unsigned int rowToRemove);//remove a row from a matrix
void removeColumn(MatrixXd& matrix, unsigned int colToRemove);//remove a column from a matrix

void  Normalization_Condition(MatrixXd &xy_i1, MatrixXd &xy_i2, MatrixXd& H1, MatrixXd& H2);

void Perform_LinOri(CameraParam& camera_params, MatrixXd &xy_i1, MatrixXd &xy_i2, MatrixXd& Fun, MatrixXd& Emat);

bool Decompose_Essential(CameraParam& camera_params, MatrixXd Emat, MatrixXd xy_i1, MatrixXd xy_i2,	double& Bx, double& By, double& Bz, double& Omega, double& Phi, double& Kappa);

bool Perform_NonlinOri(CameraParam& camera_params, MatrixXd &xy_i1, MatrixXd &xy_i2, MatrixXd& T21, double& Bx, double& By, double& Bz, double& Omega, double& Phi, double& Kappa);

void Vanilla_RANSAC(CameraParam camera_params, MatrixXd &xy_i1, MatrixXd &xy_i2, double confidence_level, double outlier_percentage, int min_iterations, double dThreshold, MatrixXd& Inliers_index);

void FindInliers(MatrixXd inliers, MatrixXd all_ties_1, MatrixXd all_ties_2, MatrixXd &inlier_ties1, MatrixXd &inlier_ties2);

MatrixXd intersection(MatrixXd x_obs_1, MatrixXd x_obs_2, RelativeOrientation ROP_1, RelativeOrientation ROP_2, CameraParam cam_params, string outfile_name);

VectorXd calculatePlane(MatrixXd Points3);

MatrixXd Compute_b(MatrixXd x_obs, MatrixXd ML, MatrixXd MR, double c, MatrixXd x_c);

MatrixXd Compute_MR(double w, double phi, double k);

MatrixXd Compute_Est(MatrixXd x_obs, MatrixXd ML, MatrixXd MR, double c);

MatrixXd Compute_A(MatrixXd x_unk, MatrixXd ML, MatrixXd MR, double c, MatrixXd x_c);

MatrixXd Compute_w(MatrixXd x_unk, MatrixXd ML, MatrixXd MR, double c, MatrixXd x_c);

MatrixXd merge_Xobs(MatrixXd x1, MatrixXd x2, CameraParam cam_params);

MatrixXd Compute_A_int(MatrixXd x_est, CameraParam params, RelativeOrientation RO1, RelativeOrientation RO2);

MatrixXd ComputeIntersectionEstimation(MatrixXd xy1, MatrixXd xy2, Matrix3d m_rot, double c, RelativeOrientation RO1, RelativeOrientation RO2, RelativeOrientation RO3);


#endif
