
#include "Lab4_Edmond.h"

double ransac_confidence = 1;
double outlier_percentage = 0.1;
int min_iterations = 2;
double dThreshold = 1;
vector<RelativeOrientation> ROs;
CameraParam camera_params;
MatrixXd inlier_temp;
MatrixXd Essential_mat;

int main() {

	//Constants
	char infile[256] = ".\\AllTies_sparsesift.txt";
	MatrixXd tie_pts, obs_01_img0, obs_02_img0, obs_01_img1, obs_12_img1, obs_02_img2, obs_12_img2;
	

	//Read in matrix
	Read_Mat(infile, tie_pts);
	// Define Camera Parameters
	camera_params.PS = 0.008609300;
	camera_params.Cn = 2592.000000000;
	camera_params.Rn = 1728.000000000;
	camera_params.xpp = 0.056678826;
	camera_params.ypp = -0.151662144;
	camera_params.f_l = 18.453986620;
	camera_params.K1 = 0.000531409551336;
	camera_params.K2 = -0.000001076228753;
	camera_params.K3 = -0.000000001031880;
	camera_params.P1 = -0.000043245177660;
	camera_params.P2 = 0.000046157068436;
	camera_params.S1 = -0.002768610675315;
	camera_params.S2 = -0.001364918105319;


	//Perform ransac on image pairs:
	SplitObs_and_RANSAC(camera_params, tie_pts, ransac_confidence, outlier_percentage, min_iterations, dThreshold);

	//Read in the ransacked image obs
	Read_Mat(".\\Inliers_Img01_0.txt", obs_01_img0);
	Read_Mat(".\\Inliers_Img02_0.txt", obs_02_img0);
	Read_Mat(".\\Inliers_Img01_1.txt", obs_01_img1);
	Read_Mat(".\\Inliers_Img12_1.txt", obs_01_img1);
	Read_Mat(".\\Inliers_Img02_2.txt", obs_02_img2);
	Read_Mat(".\\Inliers_Img12_2.txt", obs_12_img2);

	//Strip the indeces off
	MatrixXd xy_01_0 = obs_01_img0;
	MatrixXd xy_01_1 = obs_01_img0;
	MatrixXd xy_12_1 = obs_01_img0;
	MatrixXd xy_12_2 = obs_01_img0;
	MatrixXd xy_02_0 = obs_01_img0;
	MatrixXd xy_02_2 = obs_01_img0;
	removeColumn(xy_01_0,0);
	removeColumn(xy_01_1,0);
	removeColumn(xy_12_1,0);
	removeColumn(xy_12_2,0);
	removeColumn(xy_02_0,0);
	removeColumn(xy_02_2,0);
	
	//Find EOP
	MatrixXd F1, F2, F3, E1, E2, E3;
	RelativeOrientation RO1, RO2, RO3;
	RO1.bx = 0;
	RO1.by = 0;
	RO1.bz = 0; 
	RO1.kappa = 0;
	RO1.omega = 0;
	RO1.phi = 0;
	ROs.push_back(RO1);
	//First Pair
	Perform_LinOri(camera_params, xy_01_0, xy_01_1, F1, E1);
	Decompose_Essential(camera_params, E1, xy_01_0, xy_01_1, RO2.bx, RO2.by, RO2.bz, RO2.omega, RO2.phi, RO2.kappa);
	ROs.push_back(RO2);
	//Second Pair
	Perform_LinOri(camera_params, xy_02_0, xy_02_2, F2, E2);
	Decompose_Essential(camera_params, E2, xy_01_0, xy_01_1, RO3.bx, RO3.by, RO3.bz, RO3.omega, RO3.phi, RO3.kappa);
	ROs.push_back(RO3);

	//Intersection
	intersection(camera_params, xy_01_0, xy_01_1, xy_02_2,RO1, RO2, RO3);

	return 0;


}
