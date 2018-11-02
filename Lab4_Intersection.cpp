
#include "Lab4_Intersection.h"

int main() {

	//Constants
	char infile[256] = "C:\\Users\\erleahy\\Documents\\Edmond\\531\\Lab4\\Edmond_Lab4_code\\ENGO531_Lab4_Intersection\\ENGO531_Lab4_Intersection\\AllTies_sparsesift.txt\\";
	MatrixXd tie_pts;
	double ransac_confidence = 1;
	double outlier_percentage = 0.1;
	int min_iterations = 100;
	double dThreshold = 1;




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


	//Perform RANSAC
	MatrixXd inliers, inlier_temp, temp;
	MatrixXd xy_1, xy_2;


	//First indeces
	double img_number1 = tie_pts(0, 0);
	double img_number2 = tie_pts(1, 0);

	double counter = 0;

	for (int i = 0; i < tie_pts.rows()/2; i++) {

		if (tie_pts(i * 2, 1) == img_number1 && tie_pts(i * 2 + 1, 1) == img_number2) {

			if (tie_pts(i * 2, 1) == tie_pts(i * 2 + 1, 1)) {
				//IMG1 matrix 
				xy_1(counter, 0) = tie_pts(i * 2, 2);
				xy_1(counter, 1) = tie_pts(i * 2, 3);
				//IMG2 matrix 
				xy_2(counter, 0) = tie_pts(i * 2+1, 2);
				xy_2(counter, 1) = tie_pts(i * 2+1, 3);
			}

		}
		else {

			//perform RANSAC
			Vanilla_RANSAC(camera_params, xy_1, xy_2, ransac_confidence, outlier_percentage, min_iterations, dThreshold, inlier_temp);
			//Append
			temp.resize(inliers.rows() + inlier_temp.rows(), inliers.cols());
			temp << inliers, inlier_temp;
			inliers = temp;
			//Clear
			xy_1.resize()


			img_number1 = tie_pts(i*2, 0);
			img_number2 = tie_pts(i*2+1, 0);


		}
	}
	














	return 0;

}
