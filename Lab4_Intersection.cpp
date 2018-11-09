
#include "Lab4_Intersection.h"

int main() {

	//Constants
	char infile[256] = ".\\AllTies_sparsesift.txt";
	MatrixXd tie_pts, obs_01_img0, obs_02_img0, obs_01_img1, obs_12_img1, obs_02_img2, obs_12_img2;
	double ransac_confidence = 1;
	double outlier_percentage = 0.1;
	int min_iterations = 5;
	double dThreshold = 1;
	vector<int> indeces;



	//Read in matrix
	Read_Mat(infile, tie_pts);
	// Define Camera Parameters
	CameraParam camera_params;
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

	//Find EOP

	
	return 0;

}
