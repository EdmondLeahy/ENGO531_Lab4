
#include "Lab4_Edmond.h"

double ransac_confidence = 0.99;
double outlier_percentage = 0.6;
int min_iterations = 2000;
double dThreshold = 1 * pow(10,-5);
double num_pairs = 0;
vector<RelativeOrientation> ROs;
CameraParam camera_params;
MatrixXd inlier_temp;
MatrixXd Essential_mat;

int main() {

	//Constants
	char infile[256] = ".\\AllTies_sparsesift.txt";
	MatrixXd tie_pts, obs_01_img0, obs_02_img0, obs_01_img1, obs_12_img1, obs_02_img2, obs_12_img2, img_pair_indeces;
	vector<MatrixXd> all_split_obs;

	//Read in matrix
	Read_Mat(infile, tie_pts);
	//Define Camera Parameters (JEFF)
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
	camera_params.sigma_obs = 0.8;


	//----------------------Perform ransac on image pairs:

	//SplitObs(camera_params, tie_pts, all_split_obs, img_pair_indeces);

	//Ransac_All_obs(camera_params, all_split_obs, ransac_confidence, outlier_percentage, min_iterations, dThreshold);

	//--------------------- Write the matrices to files:
	string outfilename;
	char *outfilechar;

	//for (int i = 0; i < all_split_obs.size() / 2; i++) {

	//	outfilename = "Inliers_Pair" + to_string(i) + "_" + to_string(0) + ".txt";
	//	outfilechar = new char[outfilename.length() + 1];
	//	strcpy(outfilechar, outfilename.c_str());
	//	Write_Mat(outfilechar, all_split_obs[i*2], 4);

	//	outfilename = "Inliers_Pair" + to_string(i) + "_" + to_string(1) + ".txt";
	//	outfilechar = new char[outfilename.length() + 1];
	//	strcpy(outfilechar, outfilename.c_str());
	//	Write_Mat(outfilechar, all_split_obs[i * 2 + 1], 4);

	//}
	//
	// ------------------------- DEBUG ------------------------------
	MatrixXd temp_obs;
	Read_Mat("Img_pair_indices.txt", img_pair_indeces);
	Read_Mat("Inliers_Pair0_0.txt", temp_obs);
	all_split_obs.push_back(temp_obs);
	Read_Mat("Inliers_Pair0_1.txt", temp_obs);
	all_split_obs.push_back(temp_obs);
	Read_Mat("Inliers_Pair1_0.txt", temp_obs);
	all_split_obs.push_back(temp_obs);
	Read_Mat("Inliers_Pair1_1.txt", temp_obs);
	all_split_obs.push_back(temp_obs);
	Read_Mat("Inliers_Pair2_0.txt", temp_obs);
	all_split_obs.push_back(temp_obs);
	Read_Mat("Inliers_Pair2_1.txt", temp_obs);
	all_split_obs.push_back(temp_obs);
	// -------------------------------------------------------------
	
	//--------------------- Find Relative Orientation:
	vector<RelativeOrientation> ROs;
	vector<MatrixXd> F_matrices, E_matrices;
	RelativeOrientation RO_temp, Zero;
	Zero.bx = 0;
	Zero.by = 0;
	Zero.bz = 0;
	Zero.kappa = 0;
	Zero.omega = 0;
	Zero.phi = 0;
	ROs.push_back(Zero);
	MatrixXd F, E, T21;
	for (int j = 0; j < all_split_obs.size() / 2; j++) {
		// for each pair

		T21 = MatrixXd::Zero(3, 1);
		T21 << RO_temp.bx,
			 RO_temp.by,
			RO_temp.bz;
		
		Perform_LinOri(camera_params, all_split_obs[j*2], all_split_obs[j * 2 + 1], F, E);
		Decompose_Essential(camera_params, E, all_split_obs[j * 2], all_split_obs[j * 2 + 1], RO_temp.bx, RO_temp.by, RO_temp.bz, RO_temp.omega, RO_temp.phi, RO_temp.kappa);
		Perform_NonlinOri(camera_params, all_split_obs[j * 2], all_split_obs[j * 2 + 1], T21, RO_temp.bx, RO_temp.by, RO_temp.bz, RO_temp.omega, RO_temp.phi, RO_temp.kappa);
		
		cout << "R0:\n" << RO_temp.bx << endl << RO_temp.by << endl << RO_temp.bz << endl << RO_temp.omega << endl << RO_temp.phi << endl << endl;
		cout << "\nE:\n" << E << endl << "F:\n" << F << endl;
		ROs.push_back(RO_temp);
		F_matrices.push_back(F);
		E_matrices.push_back(E);

	}



	//--------------------- Perform Intersection:
	cout << "\n\n------------- START INTERSECTION ---------------\n\n";

	string int_outname1, int_outname2;
	char *x_obs_char;
	MatrixXd x_obs;
	for (int i = 0; i < all_split_obs.size() / 2; i++) {

		cout << "\n\nIntersection: " << i * 2 << " and " << i * 2 + 1 << endl;
		int_outname1 = "Intersection_Pair" + to_string((int)img_pair_indeces(i, 0)) + "_" + to_string((int)img_pair_indeces(i, 1)) + ".txt";
		x_obs = intersection(all_split_obs[i * 2], all_split_obs[i * 2 + 1], ROs[img_pair_indeces(i,0)], ROs[img_pair_indeces(i, 1)], camera_params, int_outname1);
		int_outname1 = "Intersection_Pair" + to_string((int)img_pair_indeces(i, 0)) + "_" + to_string((int)img_pair_indeces(i, 1)) + "_FINAL_Obs.txt";
		x_obs_char = &int_outname1[0u];
		Write_Mat(x_obs_char, x_obs, 7);
	}
	cout << endl;
	system("pause");
	return 0;


}
