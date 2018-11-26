
#include "Lab4_Edmond.h"

double ransac_confidence = 0.99;
double outlier_percentage = 0.6;
int min_iterations = 200;
double dThreshold = 1 * pow(10,-5);
double num_pairs = 0;
vector<RelativeOrientation> ROs;
CameraParam camera_params;
MatrixXd inlier_temp;
MatrixXd Essential_mat;

int main() {

	//Constants
	char infile[256] = ".\\AllTies_sparsesift_Kate.txt";
	MatrixXd tie_pts, obs_01_img0, obs_02_img0, obs_01_img1, obs_12_img1, obs_02_img2, obs_12_img2;
	

	//Read in matrix
	Read_Mat(infile, tie_pts);
	// Define Camera Parameters (JEFF)
	/*camera_params.PS = 0.008609300;
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
	camera_params.sigma_obs = 0.8;
*/
	// Define Camera Parameters (KATE)
	camera_params.PS = 0.00372;
	camera_params.Cn = 6000.0;
	camera_params.Rn = 4000.0;
	camera_params.xpp = 0.000165124806008;
	camera_params.ypp = - 0.000000237476973;
	camera_params.f_l = 0.000000000062796;
	camera_params.K1 = -0.000011690064009;
	camera_params.K2 = -0.000006395424825;
	camera_params.K3 = -0.000047133657656;
	camera_params.P1 = -0.000014234944030;
	camera_params.P2 = 0.066758636;
	camera_params.S1 = 0.162602957;
	camera_params.S2 = 24.736669359;
	camera_params.sigma_obs = 0.8;


	//Perform ransac on image pairs:
	SplitObs_and_RANSAC(camera_params, tie_pts, ransac_confidence, outlier_percentage, min_iterations, dThreshold, num_pairs);

	//Read in the ransacked image obs
	Read_Mat(".\\Inliers_Img01_0.txt", obs_01_img0);
	Read_Mat(".\\Inliers_Img02_0.txt", obs_02_img0);
	Read_Mat(".\\Inliers_Img01_1.txt", obs_01_img1);
	Read_Mat(".\\Inliers_Img12_1.txt", obs_01_img1);
	Read_Mat(".\\Inliers_Img02_2.txt", obs_02_img2);
	Read_Mat(".\\Inliers_Img12_2.txt", obs_12_img2);

	//Strip the indeces off
	MatrixXd xy_01_0 = obs_01_img0.rightCols(2);
	MatrixXd xy_01_1 = obs_01_img0.rightCols(2);
	MatrixXd xy_12_1 = obs_01_img0.rightCols(2);
	MatrixXd xy_12_2 = obs_01_img0.rightCols(2);
	MatrixXd xy_02_0 = obs_01_img0.rightCols(2);
	MatrixXd xy_02_2 = obs_01_img0.rightCols(2);
	
	//cout << endl << endl << xy_01_0 << endl;


	//Find EOP
	MatrixXd F1, F2, F3, E1, E2, E3;
	RelativeOrientation Zero, RO1, RO2, RO3;
	Zero.bx = 0;
	Zero.by = 0;
	Zero.bz = 0; 
	Zero.kappa = 0;
	Zero.omega = 0;
	Zero.phi = 0;
	//First Pair
	Perform_LinOri(camera_params, xy_01_0, xy_01_1, F1, E1);	
	cout << "\n\nE:\n" << E1 << "\n\nF:\n" << F1 << endl;
	Decompose_Essential(camera_params, E1, xy_01_0, xy_01_1, RO1.bx, RO1.by, RO1.bz, RO1.omega, RO1.phi, RO1.kappa);
	ROs.push_back(RO1);
	//Second Pair
	Perform_LinOri(camera_params, xy_02_0, xy_02_2, F2, E2);
	Decompose_Essential(camera_params, E2, xy_02_0, xy_02_2, RO2.bx, RO2.by, RO2.bz, RO2.omega, RO2.phi, RO2.kappa);
	ROs.push_back(RO2);
	////Third Pair
	//Perform_LinOri(camera_params, xy_12_1, xy_12_2, F3, E3);
	//Decompose_Essential(camera_params, E2, xy_12_1, xy_12_2, RO3.bx, RO3.by, RO3.bz, RO3.omega, RO3.phi, RO3.kappa);
	//ROs.push_back(RO3);

	cout << "\n\n------------- START INTERSECTION ---------------\n\n";

	cout << "R01:\n" << RO1.bx << endl << RO1.by << endl << RO1.bz << endl << RO1.omega << endl << RO1.phi << endl << endl;
	cout << "R02:\n" << RO2.bx << endl << RO2.by << endl << RO2.bz << endl << RO2.omega << endl << RO2.phi << endl << endl;

	//Intersection img 1-2
	intersection(xy_01_0, xy_01_1, Zero, RO1, camera_params, "Intersection_01.txt");


	//Intersection img 1-3
	intersection(xy_02_0, xy_02_2, Zero, RO2, camera_params, "Intersection_02.txt");


	//Intersection img 2-3
	intersection(xy_12_1, xy_12_2, RO1, RO2, camera_params, "Intersection_12.txt");


	system("pause");
	return 0;


}
