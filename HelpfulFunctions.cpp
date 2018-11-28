//In your codes you may need any or all of these libraries to be included:
/*#include <iostream>
#include <Eigen\Dense> 
#include <Eigen\Core>
#include <fstream> 
#include <stdlib.h>
#include <vector>
#include <math.h>  
#include <algorithm>
#include <stdexcept>
#include <iomanip> 
#include <string>*/

#include "Lab4_Edmond.h"

// using namespace std;

/////////////////////////////////////////
////////////////////////////////////////
//Receives the name of a file "FaileName" containing a numerical matrix, and read the matrix data into variable "m"
void Read_Mat(char *FileName, MatrixXd& m) {

	m.resize(0, 0);

	ifstream matfile;
	matfile.open(FileName, ios::in); //to open the file and start reading from the beginning
	if (matfile.fail()) //check if the file is opened successfully
	{
		cout << "\nThere was a problem reading the following file: " << endl << FileName << endl;
		//exit (EXIT_FAILURE);
		return;
	}
	else {
		cout << "\nFile read correctly. Continuing.\n";
	}

	char* readlinechr = new char[MaxMatSize];
	vector<double> v_all;
	int nrow = 0;

	while (matfile.getline(readlinechr, MaxMatSize, '\n')) {
		nrow++;
		int stln = strlen(readlinechr);
		char* readlinestr = new char[stln + 1];
		for (int i = 0; i<stln; i++)
		{
			readlinestr[i] = readlinechr[i];
		}

		readlinestr[stln] = '\0';

		stringstream rowstream(readlinestr);
		double value;
		while (!rowstream.eof()) {
			rowstream >> value;
			v_all.push_back(value);
		}
	}
	matfile.close();

	int ncol = v_all.size() / nrow;
	m.resize(nrow, ncol);

	for (int i = 0; i<nrow; i++) {
		for (int j = 0; j<ncol; j++) {
			m(i, j) = v_all.at(i*ncol + j);
		}
	}

	return;
};

//////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//Recives the name of a file "FileName", and writes the data in the matrix "m" to this file, with fixed precision "decimal_precision"
void Write_Mat(char *FileName, MatrixXd &m, int decimal_precision) {

	ofstream matfile;
	matfile.open(FileName, ios::out); //to open the file and start writing from the beginning (Notice! over writing)
	if (matfile.fail()) //check if the file is opened successfully
	{
		cout << "There was a problem reading the follwoing file: " << endl << FileName << endl;
		//exit (EXIT_FAILURE);
		return;
	}
	matfile.flags(ios::fixed);
	matfile.precision(decimal_precision);
	matfile << m;
	matfile.close();
	return;
};


//////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//Receives three rotation angles (Omega,Phi,Kappa) and returns the rotation matrix from the object to the image space (Rot_g2i)
//Note that the rotation from image to object space is the transpose of "Rot_g2i"
void Rotation_g2i(double Omega, double Phi, double Kappa, Matrix3b3 & Rot_g2i) {
	Matrix3b3 Mw0;
	Matrix3b3 Mf0;
	Matrix3b3 Mk0;

	// compute R_g_to_i and return it to the Rot_g2i

	Mw0 << 1, 0, 0,
		0, cos(Omega), sin(Omega),
		0, -sin(Omega), cos(Omega);

	Mf0 << cos(Phi), 0, -sin(Phi),
		0, 1, 0,
		sin(Phi), 0, cos(Phi);

	Mk0 << cos(Kappa), sin(Kappa), 0,
		-sin(Kappa), cos(Kappa), 0,
		0, 0, 1;

	Rot_g2i = Mk0 * Mf0*Mw0;
};

//////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//Remove row number "rowToRemove" form the matrix "Matrix".
//so, the output "matrix" has one row less than the input.
//Row indices start from 0.
void removeRow(MatrixXd& matrix, unsigned int rowToRemove) {
	unsigned int numRows = matrix.rows() - 1;
	unsigned int numCols = matrix.cols();

	if (rowToRemove < numRows) {
		matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) = matrix.block(rowToRemove + 1, 0, numRows - rowToRemove, numCols);
		matrix.conservativeResize(numRows, numCols);
	}
};


//////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//Remove column number "colToRemove" form the matrix "Matrix".
//so, the output "matrix" has one column less than the input.
//Column indices start from 0.
void removeColumn(MatrixXd& matrix, unsigned int colToRemove)
{
	unsigned int numRows = matrix.rows();
	unsigned int numCols = matrix.cols() - 1;

	if (colToRemove < numCols) {
		matrix.block(0, colToRemove, numRows, numCols - colToRemove) = matrix.block(0, colToRemove + 1, numRows, numCols - colToRemove);
		matrix.conservativeResize(numRows, numCols);
	}
};

//////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//Checked whether value "k" is included in vector "V"; if so, returns true and the index "I" at which "k" is found.
bool does_exist(vector<double> V, double k, int & I) {
	bool res = 0;
	I = -1;
	for (int i = 0; i<V.size(); i++) {
		if (k == V.at(i)) {
			I = i;
			res = 1;
			break;
		}
	}

	return res;
};

//////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//Checked whether integer value "k" is included in vector "V"; if so, returns true and the index "I" at which "k" is found.
bool does_exist_eigen(VectorXi V, int k, int & I) {
	bool res = 0;
	I = -1;
	for (int i = 0; i<V.size(); i++) {
		if (k == V(i)) {
			I = i;
			res = 1;
			break;
		}
	}

	return res;

};

//////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
// make a random permutation of the numbers in range_vec and returns it by reference
void randperm(vector<int> &range_vec)
{
	int i, j, t;
	int n = range_vec.size();
	for (i = 0; i<n; i++) {
		j = rand() % (n - i) + i;
		t = range_vec.at(j);
		range_vec.at(j) = range_vec.at(i);
		range_vec.at(i) = t;
	}
};
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//Given the rotation matrix form image to object space, "R", calculates the rotation angles (Omega,Phi,Kappa)
void Convert_R_to_Angles(Matrix3b3 R, double& Omega, double& Phi, double& Kappa) {

	double A11 = R(0, 0);
	double A12 = R(0, 1);
	double A13 = R(0, 2);
	double A21 = R(1, 0);
	double A22 = R(1, 1);
	double A23 = R(1, 2);
	double A33 = R(2, 2);

	if (A13 != 1 && A13 != -1) {

		Phi = asin(A13);
		Omega = atan2(-A23 / cos(Phi), A33 / cos(Phi));
		Kappa = atan2(-A12 / cos(Phi), A11 / cos(Phi));

	}

	else if (A13 == 1) {
		Phi = pi / 2;
		Omega = 0; //arbitrary
		Kappa = -Omega + atan2(A21, A22);

	}

	else if (A13 == -1) {
		Phi = -pi / 2;
		Omega = 0;
		Kappa = Omega + atan2(A21, A22);
	}

	return;
};
///////////////////////////////////////////////
//////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//Given corresponding points:
//       xy_i1: x in pixels, y in pixels (image observations of points on image 1)
//		 xy_i2: x in pixels, y in pixels (image observations of the same points on image 2)
			// So, the size of "xy_i1" and "xy_i2" must be the same
//Finds the conditioning homographies, H1 and H2
void  Normalization_Condition(MatrixXd &xy_i1, MatrixXd &xy_i2, MatrixXd& H1, MatrixXd& H2) {

	int n_p = xy_i1.rows();

	double cx1 = (xy_i1.block(0, 0, n_p, 1)).array().mean();
	double cy1 = (xy_i1.block(0, 1, n_p, 1)).array().mean();
	double cx2 = (xy_i2.block(0, 0, n_p, 1)).array().mean();
	double cy2 = (xy_i2.block(0, 1, n_p, 1)).array().mean();

	MatrixXd dx1 = xy_i1.block(0, 0, n_p, 1) - cx1 * (MatrixXd::Ones(n_p, 1));
	MatrixXd dy1 = xy_i1.block(0, 1, n_p, 1) - cy1 * (MatrixXd::Ones(n_p, 1));
	MatrixXd dx2 = xy_i2.block(0, 0, n_p, 1) - cx2 * (MatrixXd::Ones(n_p, 1));
	MatrixXd dy2 = xy_i2.block(0, 1, n_p, 1) - cy2 * (MatrixXd::Ones(n_p, 1));

	MatrixXd Temp = (dx1.array().square() + dy1.array().square()).array().sqrt();
	double d1 = Temp.array().mean();

	Temp = (dx2.array().square() + dy2.array().square()).array().sqrt();
	double d2 = Temp.array().mean();

	H1.setZero(3, 3);
	H2.setZero(3, 3);

	H1 << sqrt(2.0) / d1, 0, -(sqrt(2.0) / d1 * cx1),
		0, sqrt(2.0) / d1, -(sqrt(2.0) / d1 * cy1),
		0, 0, 1;
	H2 << sqrt(2.0) / d2, 0, -(sqrt(2.0) / d2 * cx2),
		0, sqrt(2.0) / d2, -(sqrt(2.0) / d2 * cy2),
		0, 0, 1;


	return;
};


//////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//Given: a structure object of type CameraParam, which contain the IOPs of the camera.  CameraP is defined as follows:
/*struct CamerParam {
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

};*/
//Given: corresponding image points (in pixels) from image one (xy_i1) and image 2 (xy_i2)
//Calculates the fundamental matrix (Fun) and the essential matrix (Emat)

void Perform_LinOri(CameraParam& camera_params, MatrixXd &xy_i1, MatrixXd &xy_i2, MatrixXd& Fun, MatrixXd& Emat) {

	//This is a scale factor for scaling the points coordinates to gain better condition number from the design matrix or the normal equations matrix
	double IO_scalfact = 0.1;
	double PS = camera_params.PS*IO_scalfact;
	double f_l = camera_params.f_l*IO_scalfact;
	double xpp = camera_params.xpp*IO_scalfact;
	double ypp = camera_params.ypp*IO_scalfact;
	double K1 = camera_params.K1 / pow(IO_scalfact, (int)2);
	double K2 = camera_params.K2 / pow(IO_scalfact, (int)4);
	double K3 = camera_params.K3 / pow(IO_scalfact, (int)6);
	double P1 = camera_params.P1 / IO_scalfact;
	double P2 = camera_params.P2 / IO_scalfact;
	double S1 = camera_params.S1;
	double S2 = camera_params.S2;
	double Cn = camera_params.Cn;
	double Rn = camera_params.Rn;

	double Cx = (-Cn / 2 + 0.5)*PS - xpp;
	double Cy = -(-Rn / 2 + 0.5)*PS - ypp;


	MatrixXd A1;
	A1.setZero(3, 3);

	MatrixXd A2;
	A2.setZero(3, 3);

	A1 << 1 + S1, S2, 0,
		0, 1, 0,
		0, 0, 1;

	A2 << PS, 0, Cx,
		0, -PS, Cy,
		0, 0, -f_l;

	MatrixXd ACalmat1;
	MatrixXd KCalmat1;

	ACalmat1.setZero(3, 3);
	ACalmat1 = A1 * A2; //This is the inverse of "K1", the calibration matrix of image 1
	double ss;
	ss = ACalmat1(2, 2);
	ACalmat1 = ACalmat1 / ss;
	KCalmat1 = ACalmat1.inverse();



	MatrixXd ACalmat2;
	MatrixXd KCalmat2;
	ACalmat2.setZero(3, 3);
	ACalmat2 = A1 * A2;//This is the inverse of "K2", calibration matrix of image 2
	ss = ACalmat2(2, 2);
	ACalmat2 = ACalmat2 / ss;
	KCalmat2 = ACalmat2.inverse();



	//find out the transformations H1 and H2 for conditioning the image coordinates
	MatrixXd H1;
	MatrixXd H2;
	Normalization_Condition(xy_i1, xy_i2, H1, H2);

	/*cout<<"H1="<<endl;
	cout<<H1<<endl;

	cout<<"H2="<<endl;
	cout<<H2<<endl;*/

	//set the coefficient matrix "A_mat"
	int n_p = xy_i1.rows();

	MatrixXd A_mat;
	A_mat.setZero(n_p, 9);

	for (int i = 0; i<n_p; i++) {
		double x1p, y1p, x2p, y2p;
		MatrixXd tmp1;
		tmp1.setZero(3, 1);
		tmp1 << xy_i1(i, 0),
			xy_i1(i, 1),
			1.0;

		MatrixXd tmp = H1 * tmp1;
		x1p = tmp(0, 0);
		y1p = tmp(1, 0);

		tmp1.setZero(3, 1);
		tmp1 << xy_i2(i, 0),
			xy_i2(i, 1),
			1.0;

		tmp = H2 * tmp1;
		x2p = tmp(0, 0);
		y2p = tmp(1, 0);

		MatrixXd a;
		a.setZero(1, 9);

		a << x1p * x2p, x1p*y2p, x1p, y1p*x2p, y2p*y1p, y1p, x2p, y2p, 1;
		A_mat.row(i) = a;
	}

	JacobiSVD<MatrixXd> svd(A_mat, ComputeFullV);
	MatrixXd V = svd.matrixV();


	Matrix<double, 9, 1> Fv = V.col(8);
	Matrix<double, 3, 3> Fcap_tmp;
	Fcap_tmp << Fv(0, 0), Fv(1, 0), Fv(2, 0),
		Fv(3, 0), Fv(4, 0), Fv(5, 0),
		Fv(6, 0), Fv(7, 0), Fv(8, 0);

	Matrix<double, 3, 3> Fcap;
	Fcap = (H1.transpose())*Fcap_tmp*H2;
	double scale = Fcap(2, 2);
	Fcap = (Fcap / scale);

	//uptohere
	JacobiSVD<MatrixXd> svd2(Fcap, ComputeFullV | ComputeFullU);
	MatrixXd V2 = svd2.matrixV();
	MatrixXd U2 = svd2.matrixU();
	MatrixXd D2 = svd2.singularValues();

	//MatrixXd temppp=U2*(D2.asDiagonal())*(V2.transpose());


	MatrixXd D_s;
	D_s.setZero(3, 1);
	D_s << D2(0, 0),
		D2(1, 0),
		0.0;
	Fun = U2 * (D_s.asDiagonal())*(V2.transpose());

	//cout << KCalmat1 << endl;

	Emat = (KCalmat1.transpose())*Fun*KCalmat2;


	return;
};
//////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//Given: a structure object of type CameraParam, which contain the IOPs of the camera.  CameraP is defined as follows:
/*struct CamerParam {
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


void SplitObs_and_RANSAC(CameraParam camera_params, MatrixXd tie_pts, double ransac_conf, double outlier_percentage, double min_iterations, double dThreshold)
{
	MatrixXd inliers, inlier_temp, temp;
	MatrixXd xy_1, xy_2, inlier_ties1, inlier_ties2, xy1_temp, xy2_temp;


	//First Photo Pair
	int img_number1 = tie_pts(0, 0);
	int img_number2 = tie_pts(1, 0);

	//For printing
	string outfilename;
	char *outfilechar;

	int counter = 0;
	for (int i = 0; i < tie_pts.rows() / 2; i++) {

		if (tie_pts(i * 2, 0) == img_number1 && tie_pts(i * 2 + 1, 0) == img_number2) {

			if (tie_pts(i * 2, 1) == tie_pts(i * 2 + 1, 1)) {
				//IMG1 matrix 
				xy_1.conservativeResize(xy_1.rows() + 1, 3);
				xy_2.conservativeResize(xy_2.rows() + 1, 3);
				xy_1(counter, 0) = tie_pts(i * 2, 1);
				xy_1(counter, 1) = tie_pts(i * 2, 2);
				xy_1(counter, 2) = tie_pts(i * 2, 3);
				//IMG2 matrix 
				xy_2(counter, 0) = tie_pts(i * 2+1, 1);
				xy_2(counter, 1) = tie_pts(i * 2+1, 2);
				xy_2(counter, 2) = tie_pts(i * 2+1, 3);
				counter++;
				fprintf(stdout, "Working on image pair %i:%i\tPoint:\t%d\r", img_number1, img_number2, counter);
			}

		}
		else {
			fprintf(stdout, "\nFinished Pair.\n\n", img_number1, img_number2, counter);
			//perform RANSAC
			xy1_temp = xy_1;
			xy2_temp = xy_2;
			removeColumn(xy1_temp, 0);
			removeColumn(xy2_temp, 0);

			Vanilla_RANSAC(camera_params, xy1_temp, xy2_temp, ransac_conf, outlier_percentage, min_iterations, dThreshold, inlier_temp);
			//Save inliers
			FindInliers(inlier_temp, xy_1, xy_2, inlier_ties1, inlier_ties2);
			//write to file
			//img1
			outfilename = "Inliers_Img" + to_string(img_number1) + to_string(img_number2) + "_" + to_string(img_number1) + ".txt";
			outfilechar = new char[outfilename.length() + 1];
			strcpy(outfilechar, outfilename.c_str());
			Write_Mat(outfilechar, inlier_ties1, 4);
			//img2
			outfilename = "Inliers_Img" + to_string(img_number1) + to_string(img_number2) + "_" + to_string(img_number2) + ".txt";
			outfilechar = new char[outfilename.length() + 1];
			strcpy(outfilechar, outfilename.c_str());
			Write_Mat(outfilechar, inlier_ties2, 4);
			//Clear
			counter = 0;
			xy_1.resize(0, 3);
			xy_2.resize(0, 3);

			//Next Photo Pair
			img_number1 = tie_pts(i * 2, 0);
			img_number2 = tie_pts(i * 2 + 1, 0);


		}
	}


	//Do one more ransac process
	fprintf(stdout, "\nFinished Pair.\n\n", img_number1, img_number2, counter);
	//perform RANSAC
	xy1_temp = xy_1;
	xy2_temp = xy_2;
	removeColumn(xy1_temp, 0);
	removeColumn(xy2_temp, 0);
	Vanilla_RANSAC(camera_params, xy1_temp, xy2_temp, ransac_conf, outlier_percentage, min_iterations, dThreshold, inlier_temp);
	//Save inliers
	FindInliers(inlier_temp, xy_1, xy_2, inlier_ties1, inlier_ties2);
	//write to file
	//img1
	outfilename = "Inliers_Img" + to_string(img_number1) + to_string(img_number2) + "_" + to_string(img_number1) + ".txt";
	outfilechar = new char[outfilename.length() + 1];
	strcpy(outfilechar, outfilename.c_str());
	Write_Mat(outfilechar, inlier_ties1, 4);
	//img2
	outfilename = "Inliers_Img" + to_string(img_number1) + to_string(img_number2) + "_" + to_string(img_number2) + ".txt";
	outfilechar = new char[outfilename.length() + 1];
	strcpy(outfilechar, outfilename.c_str());
	Write_Mat(outfilechar, inlier_ties2, 4);



}
;


};*/




//void SplitObs_and_RANSAC(CameraParam camera_params, MatrixXd tie_pts, double ransac_conf, double outlier_percentage, double min_iterations, double dThreshold, double& pair_counter)
//{
//	MatrixXd inliers, inlier_temp, temp;
//	MatrixXd xy_1, xy_2, inlier_ties1, inlier_ties2, xy1_temp, xy2_temp;
//
//
//	//First Photo Pair
//	int img_number1 = tie_pts(0, 0);
//	int img_number2 = tie_pts(1, 0);
//
//	//For printing
//	string outfilename;
//	char *outfilechar;
//
//	pair_counter = 0;
//	int counter = 0;
//	for (int i = 0; i < tie_pts.rows() / 2; i++) {
//
//		if (tie_pts(i * 2, 0) == img_number1 && tie_pts(i * 2 + 1, 0) == img_number2) {
//
//			if (tie_pts(i * 2, 1) == tie_pts(i * 2 + 1, 1)) {
//				//IMG1 matrix 
//				xy_1.conservativeResize(xy_1.rows() + 1, 3);
//				xy_2.conservativeResize(xy_2.rows() + 1, 3);
//				xy_1(counter, 0) = tie_pts(i * 2, 1);
//				xy_1(counter, 1) = tie_pts(i * 2, 2);
//				xy_1(counter, 2) = tie_pts(i * 2, 3);
//				//IMG2 matrix 
//				xy_2(counter, 0) = tie_pts(i * 2 + 1, 1);
//				xy_2(counter, 1) = tie_pts(i * 2 + 1, 2);
//				xy_2(counter, 2) = tie_pts(i * 2 + 1, 3);
//				counter++;
//				fprintf(stdout, "Working on image pair %i:%i\tPoint:\t%d\r", img_number1, img_number2, counter);
//			}
//
//		}
//		else {
//			fprintf(stdout, "\nFinished Pair.\n\n", img_number1, img_number2, counter);
//			pair_counter++;
//			//perform RANSAC
//			xy1_temp = xy_1;
//			xy2_temp = xy_2;
//			removeColumn(xy1_temp, 0);
//			removeColumn(xy2_temp, 0);
//
//			//cout << "\nXYTemp1:\n" << xy1_temp << "\nXYTemp2:\n" << xy2_temp << endl;
//
//
//			Vanilla_RANSAC(camera_params, xy1_temp, xy2_temp, ransac_conf, outlier_percentage, min_iterations, dThreshold, inlier_temp);
//			//Save inliers
//			Write_Mat("Inliers_temp.txt", inlier_temp, 4);
//			Write_Mat("xy1_temp.txt", xy1_temp, 4);
//			FindInliers(inlier_temp, xy_1, xy_2, inlier_ties1, inlier_ties2);
//			//write to file
//			//img1
//			outfilename = "Inliers_Img" + to_string(img_number1) + to_string(img_number2) + "_" + to_string(img_number1) + ".txt";
//			outfilechar = new char[outfilename.length() + 1];
//			strcpy(outfilechar, outfilename.c_str());
//			Write_Mat(outfilechar, inlier_ties1, 4);
//			//img2
//			outfilename = "Inliers_Img" + to_string(img_number1) + to_string(img_number2) + "_" + to_string(img_number2) + ".txt";
//			outfilechar = new char[outfilename.length() + 1];
//			strcpy(outfilechar, outfilename.c_str());
//			Write_Mat(outfilechar, inlier_ties2, 4);
//			//Clear
//			counter = 0;
//			xy_1.resize(0, 3);
//			xy_2.resize(0, 3);
//
//			//Next Photo Pair
//			img_number1 = tie_pts(i * 2, 0);
//			img_number2 = tie_pts(i * 2 + 1, 0);
//
//
//		}
//	}
//
//
//	//Do one more ransac process
//	fprintf(stdout, "\nFinished Pair.\n\n", img_number1, img_number2, counter);
//	pair_counter++;
//	//perform RANSAC
//	xy1_temp = xy_1;
//	xy2_temp = xy_2;
//	removeColumn(xy1_temp, 0);
//	removeColumn(xy2_temp, 0);
//	Vanilla_RANSAC(camera_params, xy1_temp, xy2_temp, ransac_conf, outlier_percentage, min_iterations, dThreshold, inlier_temp);
//	//Save inliers
//	FindInliers(inlier_temp, xy_1, xy_2, inlier_ties1, inlier_ties2);
//	//write to file
//	//img1
//	outfilename = "Inliers_Img" + to_string(img_number1) + to_string(img_number2) + "_" + to_string(img_number1) + ".txt";
//	outfilechar = new char[outfilename.length() + 1];
//	strcpy(outfilechar, outfilename.c_str());
//	Write_Mat(outfilechar, inlier_ties1, 4);
//	//img2
//	outfilename = "Inliers_Img" + to_string(img_number1) + to_string(img_number2) + "_" + to_string(img_number2) + ".txt";
//	outfilechar = new char[outfilename.length() + 1];
//	strcpy(outfilechar, outfilename.c_str());
//	Write_Mat(outfilechar, inlier_ties2, 4);
//
//
//
//}
void SplitObs(CameraParam camera_params, MatrixXd tie_pts, vector<MatrixXd>& all_split_obs, MatrixXd& img_indeces)
{
	MatrixXd inliers, inlier_temp, temp;
	MatrixXd xy_1, xy_2, inlier_ties1, inlier_ties2, xy1_temp, xy2_temp;
	img_indeces = MatrixXd(0, 0);


	//First Photo Pair
	int img_number1 = tie_pts(0, 0);
	int img_number2 = tie_pts(1, 0);

	//For printing
	string outfilename;
	char *outfilechar;
	int num_pairs = 0;
	int counter = 0;
	for (int i = 0; i < tie_pts.rows() / 2; i++) {

		if (tie_pts(i * 2, 0) == img_number1 && tie_pts(i * 2 + 1, 0) == img_number2) {

			if (tie_pts(i * 2, 1) == tie_pts(i * 2 + 1, 1)) {
				//IMG1 matrix 
				xy_1.conservativeResize(xy_1.rows() + 1, 3);
				xy_2.conservativeResize(xy_2.rows() + 1, 3);
				xy_1(counter, 0) = tie_pts(i * 2, 1);
				xy_1(counter, 1) = tie_pts(i * 2, 2);
				xy_1(counter, 2) = tie_pts(i * 2, 3);
				//IMG2 matrix 
				xy_2(counter, 0) = tie_pts(i * 2 + 1, 1);
				xy_2(counter, 1) = tie_pts(i * 2 + 1, 2);
				xy_2(counter, 2) = tie_pts(i * 2 + 1, 3);
				counter++;
				fprintf(stdout, "Working on image pair %i:%i\tPoint:\t%d\r", img_number1, img_number2, counter);
				//cout << xy_1 << endl;
			}

		}
		else {
			fprintf(stdout, "\nFinished Pair.\n\n", img_number1, img_number2, counter);
			// Pushback vector of observations
			all_split_obs.push_back(xy_1);
			all_split_obs.push_back(xy_2);
			// add the indeces to img list
			img_indeces.conservativeResize(img_indeces.rows() + 1, 2);
			img_indeces(num_pairs, 0) = img_number1;
			img_indeces(num_pairs, 1) = img_number2;
			num_pairs++;

			//Clear and set counter back to 0
			counter = 0;
			xy_1.resize(0, 3);
			xy_2.resize(0, 3);

			//Next Photo Pair
			img_number1 = tie_pts(i * 2, 0);
			img_number2 = tie_pts(i * 2 + 1, 0);


		}

	}

	fprintf(stdout, "\nFinished Pair.\n\n", img_number1, img_number2, counter);
	// Pushback vector of last observations
	all_split_obs.push_back(xy_1);
	all_split_obs.push_back(xy_2);
	// add the indeces to img list
	img_indeces.conservativeResize(img_indeces.rows() + 1, 2);
	img_indeces(num_pairs, 0) = img_number1;
	img_indeces(num_pairs, 1) = img_number2;
	num_pairs++;
}
void Ransac_All_obs(CameraParam camera_params, vector<MatrixXd>& all_split_obs, double ransac_conf, double outlier_percentage, double min_iterations, double dThreshold)
{
	int num_pairs = all_split_obs.size() / 2;
	MatrixXd xy1_temp, xy2_temp, inliers, inliers1, inliers2;
	vector<MatrixXd> all_split_obs_inliers;
	
	for (int i = 0; i < num_pairs; i++) {
		xy1_temp = all_split_obs[i * 2].rightCols(2);
		xy2_temp = all_split_obs[i * 2 + 1].rightCols(2);
		
		// Perform RANSAC, receive index of inliers
		Vanilla_RANSAC(camera_params, xy1_temp, xy2_temp, ransac_conf, outlier_percentage, min_iterations, dThreshold, inliers);
		// Remove outliers
		FindInliers(inliers, xy1_temp, xy2_temp, inliers1, inliers2);
		// Pushback new (only inliers)
		all_split_obs_inliers.push_back(inliers1);
		all_split_obs_inliers.push_back(inliers2);

	}
	// Set the split obs to inliers only
	all_split_obs = all_split_obs_inliers;

}
;



//Given: corresponding image points (in pixels) from image one (xy_i1) and image 2 (xy_i2)
//Given: the essential matrix (Emat)
//Decomposes the essential matrix to the ROPs (Bx,By,Bz,Omega,Phi,Kappa)
	//if decomposition is successful, returns true

bool Decompose_Essential(CameraParam& camera_params, MatrixXd Emat, MatrixXd xy_i1, MatrixXd xy_i2,
	double& Bx, double& By, double& Bz, double& Omega, double& Phi, double& Kappa) {

	int n_p = xy_i1.rows();

	bool succeed_done = 0;

	//Decompose the matrix E by svd
	JacobiSVD<MatrixXd> svd2(Emat, ComputeFullV | ComputeFullU);
	MatrixXd v = svd2.matrixV();
	MatrixXd u = svd2.matrixU();


	//enforce the constraints to E
	double detu = u.determinant();
	double detv = v.determinant();

	if (detu<0) {
		u = (-1.0)*u;
	}

	if (detv<0) {
		v = (-1.0)*v;
	}

	MatrixXd D_s;
	D_s.setZero(3, 1);
	D_s << 1.0,
		1.0,
		0.0;
	Emat = u * (D_s.asDiagonal())*(v.transpose());

	MatrixXd w;
	w.setZero(3, 3);
	w << 0, 1, 0,
		-1, 0, 0,
		0, 0, 1;

	MatrixXd z;
	z.setZero(3, 3);
	z << 0, 1, 0,
		-1, 0, 0,
		0, 0, 0;

	JacobiSVD<MatrixXd> svd3(Emat, ComputeFullV | ComputeFullU);
	v = svd3.matrixV();
	u = svd3.matrixU();

	MatrixXd rot1;
	MatrixXd rot2;
	MatrixXd t1;
	MatrixXd t2;

	rot1 = u * w * (v.transpose()); //R 2 to 1
	rot2 = u * (w.transpose()) * (v.transpose());

	t1 = -u.col(2); //R_eto1*B_2,1
	t2 = u.col(2);

	//4 possible choices of the camera matrix P2 based on the 2 possible
	//choices of rot and 2 possible signs of t.

	//This is a scale factor for scaling the points coordinates to gain better condition number from the design matrix or the normal equations matrix
	double IO_scalfact = 0.1;
	double PS = camera_params.PS*IO_scalfact;
	double f_l = camera_params.f_l*IO_scalfact;
	double xpp = camera_params.xpp*IO_scalfact;
	double ypp = camera_params.ypp*IO_scalfact;
	double K1 = camera_params.K1 / pow(IO_scalfact, (int)2);
	double K2 = camera_params.K2 / pow(IO_scalfact, (int)4);
	double K3 = camera_params.K3 / pow(IO_scalfact, (int)6);
	double P1 = camera_params.P1 / IO_scalfact;
	double P2 = camera_params.P2 / IO_scalfact;
	double S1 = camera_params.S1;
	double S2 = camera_params.S2;
	double Cn = camera_params.Cn;
	double Rn = camera_params.Rn;

	double Cx = (-Cn / 2 + 0.5)*PS - xpp;
	double Cy = -(-Rn / 2 + 0.5)*PS - ypp;



	MatrixXd A1;
	A1.setZero(3, 3);

	MatrixXd A2;
	A2.setZero(3, 3);

	A1 << 1 + S1, S2, 0,
		0, 1, 0,
		0, 0, 1;

	A2 << PS, 0, Cx,
		0, -PS, Cy,
		0, 0, -f_l;

	MatrixXd ACalmat1;
	ACalmat1.setZero(3, 3);
	ACalmat1 = A1 * A2; //This is the inverse of "K1", the calibration matrix of image 1
	double scale = ACalmat1(2, 2);
	ACalmat1 = ACalmat1 / scale;

	MatrixXd ACalmat2;
	ACalmat2.setZero(3, 3);
	ACalmat2 = A1 * A2;//This is the inverse of "K2", calibration matrix of image 2
	scale = ACalmat2(2, 2);
	ACalmat2 = ACalmat2 / scale;

	MatrixXd Pmat1;
	MatrixXd Temp;
	MatrixXd Pmat2;
	double signM1 = 1.0;
	double signM2 = 1.0;

	Temp.setZero(3, 4);
	Temp << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0;

	//projection matrix of image 1
	Pmat1 = (ACalmat1.inverse())*Temp;

	double detTemp = (Pmat1.block(0, 0, 3, 3)).determinant();
	if (detTemp<0) {
		signM1 = -1.0;
	}
	else {
		signM1 = 1.0;
	}

	MatrixXd Final_rotation;
	MatrixXd Final_t;

	MatrixXd signdepth2;
	MatrixXd signdepth1;

	Matrix<double, 3, 1> tmp31;
	Matrix<double, 3, 1> tmp1;
	Matrix<double, 3, 1> tmp2;

	double min_sum = 1.0;

	//find out for which solution the majority of the points are in front of the camera.
	for (int i = 0; i<4; i++) {

		MatrixXd rot, t;
		if (i == 0) {
			rot = rot1;
			t = t1;
		}
		if (i == 1) {
			rot = rot2;
			t = t2;
		}
		if (i == 2) {
			rot = rot1;
			t = t2;
		}
		if (i == 3) {
			rot = rot2;
			t = t1;
		}

		Temp.setZero(3, 4);
		Temp.block(0, 0, 3, 3) = rot.transpose(); //R 1 to 2
		Temp.block(0, 3, 3, 1) = -(rot.transpose())*t;

		//projection matrix of image 2
		Pmat2 = (ACalmat2.inverse())*Temp;

		detTemp = (Pmat2.block(0, 0, 3, 3)).determinant();
		if (detTemp<0) {
			signM2 = -1.0;
		}
		else {
			signM2 = 1.0;
		}


		signdepth1.setZero(n_p, 1);
		signdepth2.setZero(n_p, 1);

		for (int jj = 0; jj<n_p; jj++) {

			// we have to find the 3D points in this relative sense first (model intersection)
			double x1 = xy_i1(jj, 0);
			double y1 = xy_i1(jj, 1);

			double x2 = xy_i2(jj, 0);
			double y2 = xy_i2(jj, 1);

			MatrixXd A, L;
			A.setZero(4, 3);
			L.setZero(4, 1);


			A << Pmat1(0, 0) - x1 * Pmat1(2, 0), Pmat1(0, 1) - x1 * Pmat1(2, 1), Pmat1(0, 2) - x1 * Pmat1(2, 2),
				Pmat1(1, 0) - y1 * Pmat1(2, 0), Pmat1(1, 1) - y1 * Pmat1(2, 1), Pmat1(1, 2) - y1 * Pmat1(2, 2),
				Pmat2(0, 0) - x2 * Pmat2(2, 0), Pmat2(0, 1) - x2 * Pmat2(2, 1), Pmat2(0, 2) - x2 * Pmat2(2, 2),
				Pmat2(1, 0) - y2 * Pmat2(2, 0), Pmat2(1, 1) - y2 * Pmat2(2, 1), Pmat2(1, 2) - y2 * Pmat2(2, 2);


			L << -Pmat1(0, 3) + Pmat1(2, 3)*x1,
				-Pmat1(1, 3) + Pmat1(2, 3)*y1,
				-Pmat2(0, 3) + Pmat2(2, 3)*x2,
				-Pmat2(1, 3) + Pmat2(2, 3)*y2;

			tmp31 = (A.transpose()*A).ldlt().solve(A.transpose()*L);

			MatrixXd temp;
			temp.setZero(4, 1);
			temp << tmp31(0, 0),
				tmp31(1, 0),
				tmp31(2, 0),
				1.0;


			tmp1 = Pmat1 * temp;
			tmp2 = Pmat2 * temp;


			if (tmp1(2, 0)>0) {
				signdepth1(jj, 0) = signM1;
			}
			else {
				signdepth1(jj, 0) = -signM1;
			}

			if (tmp2(2, 0)>0) {
				signdepth2(jj, 0) = signM2;
			}
			else {
				signdepth2(jj, 0) = -signM2;
			}

		}

		double sumsign1 = signdepth1.array().sum();
		double sumsign2 = signdepth2.array().sum();

		if (sumsign1>min_sum && sumsign2>min_sum) {
			succeed_done = 1;
			//cout << "\n\tFrom Decomposition: solution# "<<i <<" is correct!"<< endl;
			Final_rotation = rot; //rotation form image 2 to image 1
			Final_t = t; //direction of translation vector from image 1 to image 2
			MatrixXd very_temp;
			very_temp.setZero(2, 1);
			very_temp << sumsign1, sumsign2;
			min_sum = very_temp.array().minCoeff();
			succeed_done = 1;
		}


	}


	if (succeed_done == 1) {
		//optimize the estimated ROs using the non-linear least-squares for ROP estimation 
		Convert_R_to_Angles(Final_rotation, Omega, Phi, Kappa);
		succeed_done = Perform_NonlinOri(camera_params, xy_i1, xy_i2, Final_t, Bx, By, Bz, Omega, Phi, Kappa);
	}
	


	return succeed_done;

};
//////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//Given: a structure object of type CameraParam, which contain the IOPs of the camera.  CameraP is defined as follows:
/*struct CamerParam {
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

};*/
//Given: corresponding image points (in pixels) from image one (xy_i1) and image 2 (xy_i2)
//Given: the initial ROPs, translation vector "T21" and the relative rotation angles (Omega,Phi,Kappa)
//Refines and returns the ROPs (Bx,By,Bz,Omega,Phi,Kappa) using non-linear least squares and the general parametrization of dependant images model 
       //if least squares converges, returns true
bool Perform_NonlinOri(CameraParam& camera_params, MatrixXd &xy_i1, MatrixXd &xy_i2, MatrixXd& T21, double& Bx, double& By, double& Bz, double& Omega, double& Phi, double& Kappa) {


	bool success_done = 1;

	//normalize the base vector
	double norm_T21 = sqrt(T21(0, 0)*T21(0, 0) + T21(1, 0)*T21(1, 0) + T21(2, 0)*T21(2, 0));
	Bx = T21(0, 0) / norm_T21;
	By = T21(1, 0) / norm_T21;
	Bz = T21(2, 0) / norm_T21;

	//This is a scale factor for scaling the points coordinates to gain better condition number from the design matrix or the normal equations matrix
	double IO_scalfact = 0.1;
	double PS = camera_params.PS*IO_scalfact;
	double f_l = camera_params.f_l*IO_scalfact;
	double xpp = camera_params.xpp*IO_scalfact;
	double ypp = camera_params.ypp*IO_scalfact;
	double K1 = camera_params.K1 / pow(IO_scalfact, (int)2);
	double K2 = camera_params.K2 / pow(IO_scalfact, (int)4);
	double K3 = camera_params.K3 / pow(IO_scalfact, (int)6);
	double P1 = camera_params.P1 / IO_scalfact;
	double P2 = camera_params.P2 / IO_scalfact;
	double S1 = camera_params.S1;
	double S2 = camera_params.S2;
	double Cn = camera_params.Cn;
	double Rn = camera_params.Rn;

	double Cx = (-Cn / 2 + 0.5)*PS - xpp;
	double Cy = -(-Rn / 2 + 0.5)*PS - ypp;


	//cout << " Initial ROPs:" << Omega << endl << Phi << endl << Kappa << endl << Bx << endl << By << endl << Bz << endl;
	MatrixXd A1;
	A1.setZero(3, 3);

	MatrixXd A2;
	A2.setZero(3, 3);

	A1 << 1 + S1, S2, 0,
		0, 1, 0,
		0, 0, 1;

	A2 << PS, 0, Cx,
		0, -PS, Cy,
		0, 0, -f_l;

	MatrixXd ACalmat1;
	ACalmat1.setZero(3, 3);
	ACalmat1 = A1 * A2; //This is the inverse of "K1", the calibration matrix of image 1
	double ss = ACalmat1(2, 2);
	ACalmat1 = ACalmat1 / ss;

	MatrixXd ACalmat2;
	ACalmat2.setZero(3, 3);
	ACalmat2 = A1 * A2;//This is the inverse of "K2", calibration matrix of image 2
	ss = ACalmat2(2, 2);
	ACalmat2 = ACalmat2 / ss;

	int eq_u = 6; //number of unknown
	int n_eq = xy_i1.rows();//number of equations



	int Maxiter = 10; //conditions to stop BA
	double ErrTol = 1e-15;
	int Loopcounter = 0;
	double maxcorrection = 1e20;
	MatrixXd Sigma_xcap;


	//start of the adjustment loop
	while (Loopcounter <= Maxiter && maxcorrection>ErrTol) {

		Loopcounter = Loopcounter + 1;


		//F=augmented collinearity equations
		MatrixXd Amat;
		Amat.setZero(n_eq, eq_u); //Design matrix: rond(F)/rond(X)

		MatrixXd deltaL;
		deltaL.setZero(n_eq, 1);

		MatrixXd Hconst_mat;
		Hconst_mat.setZero(1, eq_u); //matrix of constraints (on norm of base vector to be 1)

		int eqr = -1;

		Matrix3b3 R1to2;
		Rotation_g2i(Omega, Phi, Kappa, R1to2);
		Matrix3b3 R2to1;
		R2to1 = (R1to2.transpose());


		for (int pnn = 0; pnn<n_eq; pnn++) { //point number

			MatrixXd tmp1;
			tmp1.setZero(3, 1);
			double x = xy_i1(pnn, 0);
			double y = xy_i1(pnn, 1);
			double xd = (x + 0.5 - Cn / 2)*PS;
			double yd = -(y + 0.5 - Rn / 2)*PS;

			double xdd = (xd - xpp);
			double ydd = (yd - ypp);

			double r2 = xdd * xdd + ydd * ydd;

			double x1 = xdd + xdd * (K1*r2 + K2 * r2*r2 + K3 * r2*r2*r2) + P1 * (r2 + 2 * xdd*xdd) + 2 * P2*xdd*ydd + S1 * xdd + S2 * ydd;
			double y1 = ydd + ydd * (K1*r2 + K2 * (r2*r2) + K3 * (r2*r2*r2)) + P2 * (r2 + 2 * ydd*ydd) + 2 * P1*xdd*ydd;
			double z1 = -f_l;


			x = xy_i2(pnn, 0);
			y = xy_i2(pnn, 1);
			xd = (x + 0.5 - Cn / 2)*PS;
			yd = -(y + 0.5 - Rn / 2)*PS;

			xdd = (xd - xpp);
			ydd = (yd - ypp);

			r2 = xdd * xdd + ydd * ydd;

			double x22 = xdd + xdd * (K1*r2 + K2 * r2*r2 + K3 * r2*r2*r2) + P1 * (r2 + 2 * xdd*xdd) + 2 * P2*xdd*ydd + S1 * xdd + S2 * ydd;
			double y22 = ydd + ydd * (K1*r2 + K2 * (r2*r2) + K3 * (r2*r2*r2)) + P2 * (r2 + 2 * ydd*ydd) + 2 * P1*xdd*ydd;
			double z22 = -f_l;



			MatrixXd tmp;
			tmp.setZero(3, 1);
			tmp << x22, y22, z22;

			tmp = R2to1 * tmp;

			double x2 = tmp(0, 0);
			double y2 = tmp(1, 0);
			double z2 = tmp(2, 0);


			double EQ = Bz * x1*y2 - Bz * x2*y1 - By * x1*z2 + By * x2*z1 + Bx * y1*z2 - Bx * y2*z1; //F=det([Bx By Bz; xL yL zL;x'R y'R z'R]);

			double rondEQ_Bx = y1 * z2 - y2 * z1; //derivative of F w.r.t Bx

			double rondEQ_By = x2 * z1 - x1 * z2; //derivative of F w.r.t By

			double rondEQ_Bz = x1 * y2 - x2 * y1; //derivative of F w.r.t Bz

												  //derivative of F w.r.t Omega
			double rondEQ_omega = Bx * z1*z22*cos(Phi)*cos(Omega) - Bz * x1*z22*cos(Phi)*cos(Omega) - By * x1*y22*cos(Kappa)*cos(Omega) + Bx * y1*y22*cos(Kappa)*cos(Omega) + By * x1*z22*cos(Phi)*sin(Omega) - Bx * y1*z22*cos(Phi)*sin(Omega) - By * x1*x22*cos(Omega)*sin(Kappa) + Bx * x22*y1*cos(Omega)*sin(Kappa) - Bz * x1*y22*cos(Kappa)*sin(Omega) + Bx * y22*z1*cos(Kappa)*sin(Omega) - Bz * x1*x22*sin(Kappa)*sin(Omega) + Bx * x22*z1*sin(Kappa)*sin(Omega) + Bz * x1*x22*cos(Kappa)*cos(Omega)*sin(Phi) - Bx * x22*z1*cos(Kappa)*cos(Omega)*sin(Phi) - By * x1*x22*cos(Kappa)*sin(Phi)*sin(Omega) + Bx * x22*y1*cos(Kappa)*sin(Phi)*sin(Omega) - Bz * x1*y22*cos(Omega)*sin(Phi)*sin(Kappa) + Bx * y22*z1*cos(Omega)*sin(Phi)*sin(Kappa) + By * x1*y22*sin(Phi)*sin(Kappa)*sin(Omega) - Bx * y1*y22*sin(Phi)*sin(Kappa)*sin(Omega);

			////derivative of F w.r.t Phi
			double rondEQ_phi = By * z1*z22*cos(Phi) - Bz * y1*z22*cos(Phi) + Bz * x22*y1*cos(Kappa)*sin(Phi) - By * x22*z1*cos(Kappa)*sin(Phi) + By * x1*z22*cos(Omega)*sin(Phi) - Bx * y1*z22*cos(Omega)*sin(Phi) - Bz * y1*y22*sin(Phi)*sin(Kappa) + By * y22*z1*sin(Phi)*sin(Kappa) + Bz * x1*z22*sin(Phi)*sin(Omega) - Bx * z1*z22*sin(Phi)*sin(Omega) + Bz * x1*x22*cos(Phi)*cos(Kappa)*sin(Omega) - By * x1*y22*cos(Phi)*cos(Omega)*sin(Kappa) - Bx * x22*z1*cos(Phi)*cos(Kappa)*sin(Omega) + Bx * y1*y22*cos(Phi)*cos(Omega)*sin(Kappa) - Bz * x1*y22*cos(Phi)*sin(Kappa)*sin(Omega) + Bx * y22*z1*cos(Phi)*sin(Kappa)*sin(Omega) + By * x1*x22*cos(Phi)*cos(Kappa)*cos(Omega) - Bx * x22*y1*cos(Phi)*cos(Kappa)*cos(Omega);

			//derivative of F w.r.t Kappa
			double rondEQ_kappa = Bz * y1*y22*cos(Phi)*cos(Kappa) - By * y22*z1*cos(Phi)*cos(Kappa) + Bz * x1*x22*cos(Kappa)*cos(Omega) - Bx * x22*z1*cos(Kappa)*cos(Omega) + Bz * x22*y1*cos(Phi)*sin(Kappa) - By * x22*z1*cos(Phi)*sin(Kappa) - By * x1*x22*cos(Kappa)*sin(Omega) + Bx * x22*y1*cos(Kappa)*sin(Omega) - Bz * x1*y22*cos(Omega)*sin(Kappa) + Bx * y22*z1*cos(Omega)*sin(Kappa) + By * x1*y22*sin(Kappa)*sin(Omega) - Bx * y1*y22*sin(Kappa)*sin(Omega) - By * x1*y22*cos(Kappa)*cos(Omega)*sin(Phi) + Bx * y1*y22*cos(Kappa)*cos(Omega)*sin(Phi) - By * x1*x22*cos(Omega)*sin(Phi)*sin(Kappa) + Bx * x22*y1*cos(Omega)*sin(Phi)*sin(Kappa) - Bz * x1*y22*cos(Kappa)*sin(Phi)*sin(Omega) + Bx * y22*z1*cos(Kappa)*sin(Phi)*sin(Omega) - Bz * x1*x22*sin(Phi)*sin(Kappa)*sin(Omega) + Bx * x22*z1*sin(Phi)*sin(Kappa)*sin(Omega);

			eqr = eqr + 1;
			deltaL(eqr, 0) = -EQ; //negative of misclosure vector

			MatrixXd a;
			a.setZero(1, 6);
			a << rondEQ_Bx, rondEQ_By, rondEQ_Bz, rondEQ_omega, rondEQ_phi, rondEQ_kappa;
			Amat.row(eqr) = a;
		}

		//constraint equation: (sqrt(Bx*Bx + By*By + Bz*Bz))-1=0

		//the derivative of the constraint on the norm of the base vector w.r.t [Bx, By, Bz, Omega, Phi, Kapa]
		Hconst_mat << Bx / sqrt(Bx*Bx + By * By + Bz * Bz), By / sqrt(Bx*Bx + By * By + Bz * Bz), Bz / sqrt(Bx*Bx + By * By + Bz * Bz), 0, 0, 0;

		//the negative of constraint equation misclosure
		double Hc = 1 - (sqrt(Bx*Bx + By * By + Bz * Bz));

		//Temp=augmented normal equations matrix Temp=[N H';H 0]
		MatrixXd Temp;
		Temp.setZero(7, 7);
		Temp.block(0, 0, 6, 6) = (Amat.transpose())*Amat;
		Temp.block(0, 6, 6, 1) = Hconst_mat.transpose();
		Temp.block(6, 0, 1, 6) = Hconst_mat;

		//Taking the inverse of Temp and assigning it to Temp itself
		FullPivHouseholderQR<MatrixXd> qr;
		qr = Temp.fullPivHouseholderQr();

		if ((qr.rank())<7) {
			cout << "Problem here! RANK!!!!!!!!!\n";
			success_done = 0;
			break;
		}

		Temp = qr.inverse();


		//Temp_u=augmented normal equations vector Temp_u=[-A'*w;Hc];
		MatrixXd Temp_u;
		Temp_u.setZero(7, 1);
		Temp_u.block(0, 0, 6, 1) = (Amat.transpose())*deltaL;
		Temp_u(6, 0) = Hc;


		//compute the covariance matrix of unknown parameters
		Sigma_xcap = Temp.block(0, 0, 6, 6);

		//compute the corrections of parameters
		MatrixXd delta_xcap = Temp * Temp_u;

		//update the parameters
		Bx = Bx + delta_xcap(0, 0);
		By = By + delta_xcap(1, 0);
		Bz = Bz + delta_xcap(2, 0);
		Omega = Omega + delta_xcap(3, 0);
		Phi = Phi + delta_xcap(4, 0);
		Kappa = Kappa + delta_xcap(5, 0);

		maxcorrection = (delta_xcap.block(0, 0, 6, 1)).array().abs().maxCoeff();

		//show the maximum correction to the RO parameters
		//cout << "From ROP Estimation - > max deltaXcap=" << maxcorrection << endl;
		//fprintf(stdout, "From ROP Estimation - > max dletaXcap= %f\n", maxcorrection);
	}//endof the adjustment loop

	if (maxcorrection>1) {
		cout << "Problem here! Max Correction too big";
		success_done = 0;
	}

	

	return success_done;
};
//////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//Given: a structure object of type CameraParam, which contain the IOPs of the camera.  CameraP is defined as follows:
/*struct CamerParam {
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

};*/
//Given: corresponding image points (in pixels) from image one (xy_i1) and image 2 (xy_i2)
//Given: the confidence level (e.g. 0.99)
//Given: outlier percentage (e.g. 0.5) ; this is \eps in the lecture notes
//Given: minimum number of iterations before ending RANSAC
//Given: a threshold on the residuals (dThreshold) to identify the outliers from inliers
//Returns the Inliers_index, which is of the same length as "xy_i1" and "xy_i2". An elemnt of this vector
//is 1 if the point is inlier and 0 otherwise
void Vanilla_RANSAC(CameraParam camera_params, MatrixXd &xy_i1, MatrixXd &xy_i2, double confidence_level, double outlier_percentage, int min_iterations, double dThreshold, MatrixXd& Inliers_index) {

	//This is a scale factor for scaling the points coordinates to gain better condition number from the design matrix or the normal equations matrix
	double IO_scalfact = 1;
	double PS = camera_params.PS*IO_scalfact;
	double f_l = camera_params.f_l*IO_scalfact;
	double xpp = camera_params.xpp*IO_scalfact;
	double ypp = camera_params.ypp*IO_scalfact;
	double K1 = camera_params.K1 / pow(IO_scalfact, (int)2);
	double K2 = camera_params.K2 / pow(IO_scalfact, (int)4);
	double K3 = camera_params.K3 / pow(IO_scalfact, (int)6);
	double P1 = camera_params.P1 / IO_scalfact;
	double P2 = camera_params.P2 / IO_scalfact;
	double S1 = camera_params.S1;
	double S2 = camera_params.S2;
	double Cn = camera_params.Cn;
	double Rn = camera_params.Rn;

	double Cx = (-Cn / 2 + 0.5)*PS - xpp;
	double Cy = -(-Rn / 2 + 0.5)*PS - ypp;


	MatrixXd A1;
	A1.setZero(3, 3);

	MatrixXd A2;
	A2.setZero(3, 3);

	A1 << 1 + S1, S2, 0,
		0, 1, 0,
		0, 0, 1;

	A2 << PS, 0, Cx,
		0, -PS, Cy,
		0, 0, -f_l;

	MatrixXd ACalmat1;
	ACalmat1.setZero(3, 3);
	ACalmat1 = A1 * A2; //This is the inverse of "K1", the calibration matrix of image 1

	MatrixXd ACalmat2;
	ACalmat2.setZero(3, 3);
	ACalmat2 = A1 * A2;//This is the inverse of "K2", calibration matrix of image 2


	int n_p = xy_i1.rows(); //number of corresponding points

	int max_support = 0;

	Inliers_index.setZero(n_p, 1);
	MatrixXd Residuals;
	Residuals.setZero(n_p, 1);

	MatrixXd Inliers_temp;
	Inliers_temp.setZero(n_p, 1);
	MatrixXd Residuals_temp;
	Residuals_temp.setZero(n_p, 1);

	
	//MatrixXd Essential_mat; //essential matrix
	MatrixXd Essential_mat;
	Essential_mat.setZero(3, 3);

	//FIND THE MAXIMUM NUMBER OF RANSAC ITERATIONS
	int max_iter = ceil(log(1 - confidence_level) / log(1 - pow(outlier_percentage, 8.0)));
	if (max_iter<min_iterations) {
		max_iter = min_iterations;
	}


	//FORM A VECTOR WITH THE POINTS	INDICES
	vector<int> range_vec;
	for (int i = 0; i<n_p; i++) {
		range_vec.push_back(i);
	}


	int Itr = 0;

	cout << "Ransac .";
	//start ransac iterations
	while (Itr<max_iter) {
		cout << ".";
		//randomly pick 8 points
		randperm(range_vec);
		vector<int> indx;
		for (int i = 0; i<8; i++) {
			indx.push_back(range_vec.at(i));
		}

		//determine those 8 points from both images
		MatrixXd xy_i1_temp, xy_i2_temp;
		xy_i1_temp.setZero(8, 2);
		xy_i2_temp.setZero(8, 2);

		for (int i = 0; i<8; i++) {
			xy_i1_temp.row(i) = xy_i1.row(indx.at(i));
			xy_i2_temp.row(i) = xy_i2.row(indx.at(i));
		}

		//using the 8-point algorithm find the essential matrix
		MatrixXd Fmat, Emat;
		Perform_LinOri(camera_params, xy_i1_temp, xy_i2_temp, Fmat, Emat);

		double Bx, By, Bz, Omega, Phi, Kappa;
		bool did_succeed = 0;

		//decompose the essential matrix to get the (statistically optimal) RO parameters
		did_succeed = Decompose_Essential(camera_params, Emat, xy_i1_temp, xy_i2_temp, Bx, By, Bz, Omega, Phi, Kappa);


		if (did_succeed == 1) {

			//re-calculate the Essential matrix directly from the RO parameters
			MatrixXd T21;
			T21.setZero(3, 1);
			T21 << Bx,
				By,
				Bz;
			Matrix3b3 R1to2;
			Rotation_g2i(Omega, Phi, Kappa, R1to2);
			MatrixXd R2to1 = R1to2.transpose();
			MatrixXd T21x;
			T21x.setZero(3, 3);
			T21x << 0, -T21(2, 0), T21(1, 0),
				T21(2, 0), 0, -T21(0, 0),
				-T21(1, 0), T21(0, 0), 0;

			Emat = T21x * R2to1;


			MatrixXd dife;
			dife.setZero(1, 1);

			Inliers_temp.setZero(n_p, 1);

			//for all the points, determine whether they support this essential matrix or not
			for (int i = 0; i<n_p; i++) {

				MatrixXd xc_i1;
				xc_i1.setZero(3, 1);
				xc_i1 << xy_i1(i, 0),
					xy_i1(i, 1),
					1.0;
				xc_i1 = ACalmat1 * xc_i1;//calibrated image coordinates

				MatrixXd xc_i2;
				xc_i2.setZero(3, 1);
				xc_i2 << xy_i2(i, 0),
					xy_i2(i, 1),
					1.0;
				xc_i2 = ACalmat2 * xc_i2;//calibrated image coordinates


				

				//There are other sophisticated types of residuals that have the benefit of being less sensitive to noise. 
				//One of them is known as Sampson distance, which is the first-order approximation of geometric error. 
				//We want to use the Sampson distance in the implementation of RANSAC
				//COMPUTE THE SAMPSON DISTANCE AND ASSIGN IT TO VARIABLE "dd"
				
				dife = (xc_i1.transpose())*Emat*xc_i2; //xc_1'*E*xc_2

				double dd;
				dd = dife(0, 0)*dife(0, 0);

				MatrixXd L1 = Emat * xc_i2; //first epipolar line
				MatrixXd L2 = (Emat.transpose())*xc_i1; //second epipolar line

//Sampson residule
dd = dd * (1 / (L1(0, 0)*L1(0, 0) + L1(1, 0)*L1(1, 0) + L2(0, 0)*L2(0, 0) + L2(1, 0)*L2(1, 0)));

Residuals_temp(i, 0) = dd;

//check the consistency (support) for this point
if (dd <= dThreshold) {
	Inliers_temp(i, 0) = 1;
}
			}

			//support cardinality for this Essential matrix
			int support_temp = Inliers_temp.array().sum();

			//if this is the largest support so far
			if (support_temp > max_support) {
				max_support = support_temp;
				Inliers_index = Inliers_temp;// a point which is outlier gets a zero; a point which is inlier gets a one;
				Residuals = Residuals_temp;
				Essential_mat = Emat;
			}

		}
		else {
			cout << '*';
		}

		Itr = Itr + 1;

	}
	cout << endl;

	return;
};

void FindInliers(MatrixXd inliers, MatrixXd all_ties_1, MatrixXd all_ties_2, MatrixXd &inlier_ties1, MatrixXd &inlier_ties2)
{
	inlier_ties1 = all_ties_1;
	inlier_ties2 = all_ties_2;
	int num_removed = 0;
	if (inliers.rows() != all_ties_1.rows() || inliers.rows() != all_ties_2.rows()) {
		cout << "\nERROR:The rows for inlier detection are not the same!\n";
		return;

	}
	else {

		for (int i = 0; i < inliers.rows(); i++) {

			if (inliers(i, 0) != 1) { //point is an outlier
				removeRow(inlier_ties1, i-num_removed);
				removeRow(inlier_ties2, i-num_removed);
				num_removed++;
			}
		}

	}
};


MatrixXd intersection(MatrixXd x_obs_1, MatrixXd x_obs_2, RelativeOrientation ROP_1, RelativeOrientation ROP_2, CameraParam cam_params, string outfile_name) {
	int n_points = x_obs_1.rows(); //number of points in the intersection


	//EOP
	MatrixXd EOP(6,2); //Both EOP in one matrix (Xc, Yc, Zc, omega, phi, kappa)
	EOP(0, 0) = ROP_1.bx; 
	EOP(1, 0) = ROP_1.by; 
	EOP(2, 0) = ROP_1.bz; 
	EOP(3, 0) = ROP_1.omega;
	EOP(4, 0) = ROP_1.phi;
	EOP(5, 0) = ROP_1.kappa;

	EOP(0, 1) = ROP_2.bx;
	EOP(1, 1) = ROP_2.by;
	EOP(2, 1) = ROP_2.bz;
	EOP(3, 1) = ROP_2.omega;
	EOP(4, 1) = ROP_2.phi;
	EOP(5, 1) = ROP_2.kappa;

	// x_obs
	//combine the observations in correct manner
	MatrixXd x_obs(2*n_points,2); 
	x_obs = merge_Xobs(x_obs_1, x_obs_2);


	//output file
	ofstream out;
	out.open(outfile_name);

	out << "------------------- ROP's -------------------\nROP:\n" << EOP << endl;

	//create the matrix for the unknowns, x,y,z for each point
	MatrixXd x_unk = MatrixXd::Zero(3, x_obs.rows() / 2);

	MatrixXd P = MatrixXd::Zero(4, 4);
	P(0, 0) = 1 / pow(cam_params.sigma_obs, 2);
	P(1, 1) = 1 / pow(cam_params.sigma_obs, 2);
	P(2, 2) = 1 / pow(cam_params.sigma_obs, 2);
	P(3, 3) = 1 / pow(cam_params.sigma_obs, 2);


	//cout << pow(cam_params.sigma_obs, 2) << "	" << cam_params.sigma_obs << "	" << endl << P << endl << endl;
	MatrixXd Iden = MatrixXd::Identity(3, 3);
	MatrixXd A = MatrixXd::Zero(4, 3);
	MatrixXd w = MatrixXd::Zero(4, 1);
	MatrixXd dx = MatrixXd::Zero(3, 1);

	//cout << x_obs << endl << endl;

	MatrixXd ComputedPoints(x_obs.rows() / 2, 3);
	double itcount = 0;
	for (int i = 0; i < x_obs.rows(); i = i + 2)//for all points
	{
		fprintf(stdout, "Intersecting point:\t%d\r", (i + 2) / 2);
		out << "------------------- " << "Point: " << (i + 2) / 2 << " -----------------" << endl << endl;
		MatrixXd x_obs_pt = MatrixXd::Zero(4, 1);
		x_obs_pt << x_obs(i, 0),
			x_obs(i, 1),
			x_obs(i + 1, 0),
			x_obs(i + 1, 1);
		//cout << x_obs_pt << endl;
		//compute the rotation matrices for the Left image and the right image for the first point
		MatrixXd ML = Compute_MR(EOP(3, 0)*pi / 180, EOP(4, 0)*pi / 180, EOP(5, 0)*pi / 180);
		MatrixXd MR = Compute_MR(EOP(3, 1)*pi / 180, EOP(4, 1)*pi / 180, EOP(5, 1)*pi / 180);
		//compute design matrix for the point
		MatrixXd Est = Compute_Est(x_obs_pt, ML, MR, cam_params.f_l);	

		//compute b matrix
		MatrixXd x_c = MatrixXd::Zero(2, 3);
		x_c << EOP(0, 0), EOP(1, 0), EOP(2, 0),
			EOP(0, 1), EOP(1, 1), EOP(2, 1);
		MatrixXd b = Compute_b(x_obs_pt, ML, MR, cam_params.f_l, x_c);

		//compute initial estimate of the X, Y and Z of the point
		x_unk = (Est.transpose()*Est).inverse()*Est.transpose()*b;
		//cout << fixed << setprecision(7) <<"Approximate Values: " <<endl<< x_unk << endl<<endl;
		out << fixed << setprecision(7) << "Approximate Values: " << endl << x_unk << endl << endl;
		dx(0, 0) = 10000;
		itcount = 0;
		while (sqrt(pow(dx(0, 0), 2))>0.0000001 || sqrt(pow(dx(1, 0), 2))>0.0000001 || sqrt(pow(dx(2, 0), 2))>0.0000001 )
		{	
			if (itcount > 100) {
				out << " \nERROR: DOES NOT CONVERTGE.\n";
				break;
			}
			itcount++;
			//cout << x_unk << endl << endl;
			A = Compute_A(x_unk, ML, MR, cam_params.f_l, x_c);
			//cout << fixed << setprecision(7) <<"Intermediate A matrix"<<endl<< A << endl << endl;
			out << fixed << setprecision(7) << "Intermediate A matrix" << endl << A << endl << endl;
			//w = x_obs_pt - A*x_unk;
			MatrixXd w = x_obs_pt - Compute_w(x_unk, ML, MR, cam_params.f_l, x_c);
			out << fixed << setprecision(7) << "Intermediate w matrix" << endl << w << endl << endl;
			dx = -1 * (A.transpose()*P*A).inverse()*A.transpose()*P*w;
			out<<fixed<<setprecision(7) << dx << endl << endl;

			x_unk = x_unk - dx;
		}

		//Compute covariance matrix
		MatrixXd Cx = (A.transpose()*P*A).inverse();

		//cout << fixed << setprecision(15) << "Cx matrix" << endl << Cx << endl << endl << "Standard deviations:" << endl << sqrt(Cx(0, 0)) << endl << sqrt(Cx(1, 1)) << endl << sqrt(Cx(2, 2)) << endl << endl;
		out << fixed << setprecision(6) << "Cx matrix" << endl << Cx << endl << endl << "Standard deviations:" << endl << sqrt(Cx(0, 0)) << endl << sqrt(Cx(1, 1)) << endl << sqrt(Cx(2, 2)) << endl << endl;

		//compute redundancy numbers
		//cout << P << endl << endl;
		MatrixXd redun_num = (A*Cx*A.transpose()*P);
		//cout << fixed << setprecision(6) << "Redundancy Numbers: " << endl << redun_num << endl << 1 - redun_num(0, 0) << endl << 1 - redun_num(1, 1) << endl << 1 - redun_num(2, 2) << endl << 1 - redun_num(3, 3) << endl;
		out << fixed << setprecision(6) << "Redundancy Numbers: " << endl << redun_num << endl << 1 - redun_num(0, 0) << endl << 1 - redun_num(1, 1) << endl << 1 - redun_num(2, 2) << endl << 1 - redun_num(3, 3) << endl << 1 - redun_num(0, 0) + 1 - redun_num(1, 1) + 1 - redun_num(2, 2) + 1 - redun_num(3, 3) << endl;

		//cout << "Final A matrix: "<<endl << A << endl << endl << "Final x_unk: " <<endl<< x_unk << endl << endl;// << "Final << x_obs_pt << endl << endl;
		out << "Final x_unk: " << endl << x_unk << endl << endl;// << "Final << x_obs_pt << endl << endl;

																//compute residuals
		w = x_obs_pt - Compute_w(x_unk, ML, MR, cam_params.f_l, x_c);
		MatrixXd v = redun_num * w - w;
		//cout << "Residuals: " << endl << v << endl << endl << "----------" << endl;
		out << "Residuals: " << endl << v << endl << endl << "---------------------------------------------" << endl;

		//Place final points in Matrix
		ComputedPoints(i / 2, 0) = x_unk(0);
		ComputedPoints(i / 2, 1) = x_unk(1);
		ComputedPoints(i / 2, 2) = x_unk(2);
	}



	//cout << "\n----------\n" << "FINAL POINTS: \n" << ComputedPoints << "\n\n----------\n";
	out << "----------\n" << "FINAL POINTS: \n" << ComputedPoints << "\n\n----------\n";

	out.close();

	return ComputedPoints;

}
VectorXd calculatePlane(MatrixXd Points3)
{
	int num_points = Points3.rows();
	//Vector of results
	VectorXd planeParams(3);
	//Design matrix
	MatrixXd A(num_points, 3);
	//Observation Matrix
	VectorXd w(num_points, 1);


	for (int i = 0; i < num_points; i++) {
		A(i, 0) = Points3(i, 0);
		A(i, 1) = Points3(i, 1);
		A(i, 2) = 1;
		w(i, 0) = Points3(i, 2);
	}
	MatrixXd N = (A.transpose()*A).inverse();
	MatrixXd U = (A.transpose()*w);

	planeParams = N * U;

	return planeParams;
}
;

MatrixXd Compute_b(MatrixXd x_obs, MatrixXd ML, MatrixXd MR, double c, MatrixXd x_c)
{
	//x_obs should be image pt of left image x, y, right image x,y
	//M should be the rotation matrix
	//c should be the persepctive distance
	//x_c should be the left image Xc, Yc, Zc; right image Xc, Yc, Zc

	MatrixXd b = MatrixXd::Zero(4, 1);

	b << (x_obs(0, 0)*ML(2, 0) + c * ML(0, 0))*x_c(0, 0) + (x_obs(0, 0)*ML(2, 1) + c * ML(0, 1))*x_c(0, 1) + (x_obs(0, 0)*ML(2, 2) + c * ML(0, 2))*x_c(0, 2),
		(x_obs(1, 0)*ML(2, 0) + c * ML(1, 0))*x_c(0, 0) + (x_obs(1, 0)*ML(2, 1) + c * ML(1, 1))*x_c(0, 1) + (x_obs(1, 0)*ML(2, 2) + c * ML(1, 2))*x_c(0, 2),
		(x_obs(2, 0)*MR(2, 0) + c * MR(0, 0))*x_c(1, 0) + (x_obs(2, 0)*MR(2, 1) + c * MR(0, 1))*x_c(1, 1) + (x_obs(2, 0)*MR(2, 2) + c * MR(0, 2))*x_c(1, 2),
		(x_obs(3, 0)*MR(2, 0) + c * MR(1, 0))*x_c(1, 0) + (x_obs(3, 0)*MR(2, 1) + c * MR(1, 1))*x_c(1, 1) + (x_obs(3, 0)*MR(2, 2) + c * MR(1, 2))*x_c(1, 2);

	return b;
}

MatrixXd Compute_MR(double w, double phi, double k)
{
	MatrixXd MR(3, 3);
	//double w = x0(0, 0);
	//double phi = x0(1, 0);
	//double k = x0(2, 0);
	MR(0, 0) = cos(phi)*cos(k);
	MR(0, 1) = cos(w)*sin(k) + sin(w)*sin(phi)*cos(k);
	MR(0, 2) = sin(w)*sin(k) - cos(w)*sin(phi)*cos(k);
	MR(1, 0) = -cos(phi)*sin(k);
	MR(1, 1) = cos(w)*cos(k) - sin(w)*sin(phi)*sin(k);
	MR(1, 2) = sin(w)*cos(k) + cos(w)*sin(phi)*sin(k);
	MR(2, 0) = sin(phi);
	MR(2, 1) = -sin(w)*cos(phi);
	MR(2, 2) = cos(w)*cos(phi);

	return MR;
}

MatrixXd Compute_Est(MatrixXd x_obs, MatrixXd ML, MatrixXd MR, double c)
{
	//x_obs should be image pts of left image x, y, right image x,y
	//M should be the rotation matrix
	//c should be the persepctive distance

	MatrixXd Est = MatrixXd::Zero(4, 3);

	Est << x_obs(0, 0)*ML(2, 0) + c * ML(0, 0), x_obs(0, 0)*ML(2, 1) + c * ML(0, 1), x_obs(0, 0)*ML(2, 2) + c * ML(0, 2),
		x_obs(1, 0)*ML(2, 0) + c * ML(1, 0), x_obs(1, 0)*ML(2, 1) + c * ML(1, 1), x_obs(1, 0)*ML(2, 2) + c * ML(1, 2),
		x_obs(2, 0)*MR(2, 0) + c * MR(0, 0), x_obs(2, 0)*MR(2, 1) + c * MR(0, 1), x_obs(2, 0)*MR(2, 2) + c * MR(0, 2),
		x_obs(3, 0)*MR(2, 0) + c * MR(1, 0), x_obs(3, 0)*MR(2, 1) + c * MR(1, 1), x_obs(3, 0)*MR(2, 2) + c * MR(1, 2);


	return Est;
}

MatrixXd Compute_A(MatrixXd x_unk, MatrixXd ML, MatrixXd MR, double c, MatrixXd x_c)
{
	MatrixXd A = MatrixXd::Zero(4, 3);

	double UL = ML(0, 0)*(x_unk(0, 0) - x_c(0, 0)) + ML(0, 1)*(x_unk(1, 0) - x_c(0, 1)) + ML(0, 2)*(x_unk(2, 0) - x_c(0, 2));
	double WL = ML(2, 0)*(x_unk(0, 0) - x_c(0, 0)) + ML(2, 1)*(x_unk(1, 0) - x_c(0, 1)) + ML(2, 2)*(x_unk(2, 0) - x_c(0, 2));
	double VL = ML(1, 0)*(x_unk(0, 0) - x_c(0, 0)) + ML(1, 1)*(x_unk(1, 0) - x_c(0, 1)) + ML(1, 2)*(x_unk(2, 0) - x_c(0, 2));

	double UR = MR(0, 0)*(x_unk(0, 0) - x_c(1, 0)) + MR(0, 1)*(x_unk(1, 0) - x_c(1, 1)) + MR(0, 2)*(x_unk(2, 0) - x_c(1, 2));
	double WR = MR(2, 0)*(x_unk(0, 0) - x_c(1, 0)) + MR(2, 1)*(x_unk(1, 0) - x_c(1, 1)) + MR(2, 2)*(x_unk(2, 0) - x_c(1, 2));
	double VR = MR(1, 0)*(x_unk(0, 0) - x_c(1, 0)) + MR(1, 1)*(x_unk(1, 0) - x_c(1, 1)) + MR(1, 2)*(x_unk(2, 0) - x_c(1, 2));

	A << (c / pow(WL, 2))*(ML(2, 0)*UL - ML(0, 0)*WL), (c / pow(WL, 2))*(ML(2, 1)*UL - ML(0, 1)*WL), (c / pow(WL, 2))*(ML(2, 2)*UL - ML(0, 2)*WL),
		(c / pow(WL, 2))*(ML(2, 0)*VL - ML(1, 0)*WL), (c / pow(WL, 2))*(ML(2, 1)*VL - ML(1, 1)*WL), (c / pow(WL, 2))*(ML(2, 2)*VL - ML(1, 2)*WL),
		(c / pow(WR, 2))*(MR(2, 0)*UR - MR(0, 0)*WR), (c / pow(WR, 2))*(MR(2, 1)*UR - MR(0, 1)*WR), (c / pow(WR, 2))*(MR(2, 2)*UR - MR(0, 2)*WR),
		(c / pow(WR, 2))*(MR(2, 0)*VR - MR(1, 0)*WR), (c / pow(WR, 2))*(MR(2, 1)*VR - MR(1, 1)*WR), (c / pow(WR, 2))*(MR(2, 2)*VR - MR(1, 2)*WR);


	return A;
}

MatrixXd Compute_w(MatrixXd x_unk, MatrixXd ML, MatrixXd MR, double c, MatrixXd x_c)
{
	MatrixXd w = MatrixXd::Zero(4, 1);

	double UL = ML(0, 0)*(x_unk(0, 0) - x_c(0, 0)) + ML(0, 1)*(x_unk(1, 0) - x_c(0, 1)) + ML(0, 2)*(x_unk(2, 0) - x_c(0, 2));
	double WL = ML(2, 0)*(x_unk(0, 0) - x_c(0, 0)) + ML(2, 1)*(x_unk(1, 0) - x_c(0, 1)) + ML(2, 2)*(x_unk(2, 0) - x_c(0, 2));
	double VL = ML(1, 0)*(x_unk(0, 0) - x_c(0, 0)) + ML(1, 1)*(x_unk(1, 0) - x_c(0, 1)) + ML(1, 2)*(x_unk(2, 0) - x_c(0, 2));

	double UR = MR(0, 0)*(x_unk(0, 0) - x_c(1, 0)) + MR(0, 1)*(x_unk(1, 0) - x_c(1, 1)) + MR(0, 2)*(x_unk(2, 0) - x_c(1, 2));
	double WR = MR(2, 0)*(x_unk(0, 0) - x_c(1, 0)) + MR(2, 1)*(x_unk(1, 0) - x_c(1, 1)) + MR(2, 2)*(x_unk(2, 0) - x_c(1, 2));
	double VR = MR(1, 0)*(x_unk(0, 0) - x_c(1, 0)) + MR(1, 1)*(x_unk(1, 0) - x_c(1, 1)) + MR(1, 2)*(x_unk(2, 0) - x_c(1, 2));

	w << -c * UL / WL,
		-c * VL / WL,
		-c * UR / WR,
		-c * VR / WR;


	return w;
}

MatrixXd merge_Xobs(MatrixXd x1, MatrixXd x2)
{
	MatrixXd x_obs(2*x1.rows(),2);

	for (int i = 0; i < x1.rows(); i++) {
		x_obs(i * 2, 0) = x1(i, 0);   //First image obs
		x_obs(i * 2, 1) = x1(i, 1);

		x_obs(i * 2 + 1, 0) = x2(i, 0); //second image obs
		x_obs(i * 2 + 1, 1) = x2(i, 1);

	}
	
	return x_obs;
}

MatrixXd Compute_A_int(MatrixXd x_est, CameraParam params, RelativeOrientation RO1, RelativeOrientation RO2) {
	MatrixXd A = MatrixXd::Zero(x_est.rows() * 2, 3);
	Matrix3d ML, MR;
	
	Rotation_g2i(RO1.omega, RO1.phi, RO1.kappa, ML);
	Rotation_g2i(RO2.omega, RO2.phi, RO2.kappa, MR);

	double c = params.Cn; //Same for both camera


	for (int i = 0; i < x_est.rows()/2; i++) {
		double UL = ML(0, 0)*(x_est(i*2, 0) - RO1.bx) + ML(0, 1)*(x_est(i*+1, 0) - RO1.by) + ML(0, 2)*(x_est(i, 0) - RO1.bz);
		double WL = ML(2, 0)*(x_est(i, 0) - RO1.bx) + ML(2, 1)*(x_est(i+1, 0) - RO1.by) + ML(2, 2)*(x_est(i, 0) - RO1.bz);
		double VL = ML(1, 0)*(x_est(i, 0) - RO1.bx) + ML(1, 1)*(x_est(i+1, 0) - RO1.by) + ML(1, 2)*(x_est(i, 0) - RO1.bz);
																	 									  
		double UR = MR(0, 0)*(x_est(i, 0) - RO2.bx) + MR(0, 1)*(x_est(i+1, 0) - RO2.by) + MR(0, 2)*(x_est(i+1, 0) - RO2.bz);
		double WR = MR(2, 0)*(x_est(i, 0) - RO2.bx) + MR(2, 1)*(x_est(i+1, 0) - RO2.by) + MR(2, 2)*(x_est(i+1, 0) - RO2.bz);
		double VR = MR(1, 0)*(x_est(i, 0) - RO2.bx) + MR(1, 1)*(x_est(i+1, 0) - RO2.by) + MR(1, 2)*(x_est(i+1, 0) - RO2.bz);

		A(i * 4, 0) = (c / pow(WL, 2))*(ML(2, 0)*UL - ML(0, 0)*WL);
		A(i * 4, 1) = (c / pow(WL, 2))*(ML(2, 1)*UL - ML(0, 1)*WL);
		A(i * 4, 2) = (c / pow(WL, 2))*(ML(2, 2)*UL - ML(0, 2)*WL);
		A(i * 4 + 1, 0) = (c / pow(WL, 2))*(ML(2, 0)*VL - ML(1, 0)*WL);
		A(i * 4 + 1, 1) = (c / pow(WL, 2))*(ML(2, 1)*VL - ML(1, 1)*WL);
		A(i * 4 + 1, 2) = (c / pow(WL, 2))*(ML(2, 2)*VL - ML(1, 2)*WL);
		A(i * 4 + 2, 0) = (c / pow(WR, 2))*(MR(2, 0)*UR - MR(0, 0)*WR);
		A(i * 4 + 2, 1) = (c / pow(WR, 2))*(MR(2, 1)*UR - MR(0, 1)*WR);
		A(i * 4 + 2, 2) = (c / pow(WR, 2))*(MR(2, 2)*UR - MR(0, 2)*WR);
		A(i * 4 + 3, 0) = (c / pow(WR, 2))*(MR(2, 0)*VR - MR(1, 0)*WR);
		A(i * 4 + 3, 1) = (c / pow(WR, 2))*(MR(2, 1)*VR - MR(1, 1)*WR);
		A(i * 4 + 3, 2) = (c / pow(WR, 2))*(MR(2, 2)*VR - MR(1, 2)*WR);

	}
	return A;

}

MatrixXd ComputeIntersectionEstimation(MatrixXd xy1, MatrixXd xy2, Matrix3d m_rot, double c, RelativeOrientation RO1, RelativeOrientation RO2, RelativeOrientation RO3) {

	MatrixXd A = MatrixXd::Zero(6 * xy1.rows(), 3);
	MatrixXd b = MatrixXd::Zero(6 * xy1.rows(), 1);
	double Xc1 = RO1.bx;
	double Yc1 = RO1.bz;
	double Zc1 = RO1.by;
	double Xc2 = RO2.bx;
	double Yc2 = RO2.bz;
	double Zc2 = RO2.by;
	double Xc3 = RO3.bx;
	double Yc3 = RO3.bz;
	double Zc3 = RO3.by;

	for (int i = 0; i < xy1.rows(); i++){

		//A matrix
		A(i * 6, 0) = (xy1(i, 0)*m_rot(2, 0) + c*m_rot(0, 0)); //a11
		A(i * 6, 1) = (xy1(i, 0)*m_rot(2, 1) + c*m_rot(0, 1));//a12
		A(i * 6, 2) = (xy1(i, 0)*m_rot(2, 2) + c*m_rot(0, 2));//a13
		A(i * 6 + 1, 0) = (xy1(i, 1)*m_rot(2, 0) + c*m_rot(1, 0)); //a21
		A(i * 6 + 1, 1) = (xy1(i, 1)*m_rot(2, 1) + c*m_rot(1, 1));//a22
		A(i * 6 + 1, 2) = (xy1(i, 1)*m_rot(2, 2) + c*m_rot(1, 2));//a23
		//Second Image
		A(i * 6 + 2, 0) = (xy2(i, 0)*m_rot(2, 0) + c*m_rot(0, 0)); //a31
		A(i * 6 + 2, 1) = (xy2(i, 0)*m_rot(2, 1) + c*m_rot(0, 1));//a32
		A(i * 6 + 2, 2) = (xy2(i, 0)*m_rot(2, 2) + c*m_rot(0, 2));//a33
		A(i * 6 + 3, 0) = (xy2(i, 1)*m_rot(2, 0) + c*m_rot(1, 0)); //a41
		A(i * 6 + 3, 1) = (xy2(i, 1)*m_rot(2, 1) + c*m_rot(1, 1));//a42
		A(i * 6 + 3, 2) = (xy2(i, 1)*m_rot(2, 2) + c*m_rot(1, 2));//a43
	    //Third Image
		A(i * 6 + 2, 0) = (xy2(i, 0)*m_rot(2, 0) + c*m_rot(0, 0)); //a31
		A(i * 6 + 2, 1) = (xy2(i, 0)*m_rot(2, 1) + c*m_rot(0, 1));//a32
		A(i * 6 + 2, 2) = (xy2(i, 0)*m_rot(2, 2) + c*m_rot(0, 2));//a33
		A(i * 6 + 3, 0) = (xy2(i, 1)*m_rot(2, 0) + c*m_rot(1, 0)); //a41
		A(i * 6 + 3, 1) = (xy2(i, 1)*m_rot(2, 1) + c*m_rot(1, 1));//a42
		A(i * 6 + 3, 2) = (xy2(i, 1)*m_rot(2, 2) + c*m_rot(1, 2));//a43
		//b Vector
		b(i * 6) = A(i * 6, 0)*Xc1 + A(i * 6, 1)*Yc1 + A(i * 6, 2)*Zc1; //b1
		b(i * 6 + 1) = A(i * 6 + 1, 0)*Xc1 + A(i * 6 + 1, 1)*Yc1 + A(i * 6 + 1, 2)*Zc1;//b2
		b(i * 6 + 2) = A(i * 6 + 2, 0)*Xc2 + A(i * 6 + 2, 1)*Yc2 + A(i * 6 + 2, 2)*Zc2; //b3
		b(i * 6 + 3) = A(i * 6 + 3, 0)*Xc2 + A(i * 6 + 3, 1)*Yc2 + A(i * 6 + 3, 2)*Zc2;//b4
		b(i * 6 + 4) = A(i * 6 + 4, 0)*Xc2 + A(i * 6 + 4, 1)*Yc2 + A(i * 6 + 4, 2)*Zc2; //b3
		b(i * 6 + 5) = A(i * 6 + 5, 0)*Xc2 + A(i * 6 + 5, 1)*Yc2 + A(i * 6 + 5, 2)*Zc2;//b4
	}

	return (A.transpose()*A).ldlt().solve(A.transpose()*b);

}

//RelativeOrientation Approx_Planar(CameraParam cam_params, MatrixXd img_points, MatrixXd obj_points)
//{
//	RelativeOrientation ROP;
//	double x, y, X, Y, Z;
//	//STEP 1: DLT Design Matrix
//	MatrixXd A = MatrixXd::Zero(2 * img_points.rows(), 9);
//	for (int i = 0; i < img_points.rows(); i++) {
//		x = img_points(i, 1);
//		y = img_points(i, 2);
//		X = obj_points(i, 1);
//		Y = obj_points(i, 2);
//		A(i * 2, 0) = X;
//		A(i * 2, 1) = Y;
//		A(i * 2, 2) = 1;
//		A(i * 2, 3) = 0;
//		A(i * 2, 4) = 0;
//		A(i * 2, 5) = 0;
//		A(i * 2, 6) = -x*X;
//		A(i * 2, 7) = -x*Y;
//		A(i * 2, 7) = -x;
//
//		A(i * 2+1, 0) = 0;
//		A(i * 2+1, 1) = 0;
//		A(i * 2+1, 2) = 0;
//		A(i * 2+1, 3) = X;
//		A(i * 2+1, 4) = Y;
//		A(i * 2+1, 5) = 1;
//		A(i * 2+1, 6) = -y*X;
//		A(i * 2+1, 7) = -y*Y;
//		A(i * 2+1, 7) = -y;
//
//	}
//
//
//	//STEP 2: Perform SV Decomp
//
//	JacobiSVD<MatrixXd> svd(A, ComputeFullV);
//	MatrixXd V = svd.matrixV();
//
//	return ROP;
//
//
//}

MatrixXd Compute_w_int(MatrixXd estimated_xy, MatrixXd xy_1_obs, MatrixXd xy_2_obs, CameraParam params, RelativeOrientation RO1, RelativeOrientation RO2)
{
	//Computes the misclosure matrix for the intersection
	// size: 4*n_points by 1
	MatrixXd w = MatrixXd::Zero(estimated_xy.rows()*4, 1);

	Matrix3d ML, MR;

	Rotation_g2i(RO1.omega, RO1.phi, RO1.kappa, ML);
	Rotation_g2i(RO2.omega, RO2.phi, RO2.kappa, MR);

	double c = params.Cn; //Same for both camera
	double point = 0;
	for (int i = 0; i < estimated_xy.rows(); i++) {
		double UL = ML(0, 0)*(estimated_xy(i, 0) - RO1.bx) + ML(0, 1)*(estimated_xy(i + 1, 0) - RO1.by) + ML(0, 2)*(estimated_xy(i, 0) - RO1.bz);
		double WL = ML(2, 0)*(estimated_xy(i, 0) - RO1.bx) + ML(2, 1)*(estimated_xy(i + 1, 0) - RO1.by) + ML(2, 2)*(estimated_xy(i, 0) - RO1.bz);
		double VL = ML(1, 0)*(estimated_xy(i, 0) - RO1.bx) + ML(1, 1)*(estimated_xy(i + 1, 0) - RO1.by) + ML(1, 2)*(estimated_xy(i, 0) - RO1.bz);

		double UR = MR(0, 0)*(estimated_xy(i, 0) - RO2.bx) + MR(0, 1)*(estimated_xy(i + 1, 0) - RO2.by) + MR(0, 2)*(estimated_xy(i + 1, 0) - RO2.bz);
		double WR = MR(2, 0)*(estimated_xy(i, 0) - RO2.bx) + MR(2, 1)*(estimated_xy(i + 1, 0) - RO2.by) + MR(2, 2)*(estimated_xy(i + 1, 0) - RO2.bz);
		double VR = MR(1, 0)*(estimated_xy(i, 0) - RO2.bx) + MR(1, 1)*(estimated_xy(i + 1, 0) - RO2.by) + MR(1, 2)*(estimated_xy(i + 1, 0) - RO2.bz);

		w(i * 4) = xy_1_obs(i*2,1) -c*UL / WL;
		w(i * 4 + 1) = -c*VL / WL;
		w(i * 4 + 2) = -c*VL / WL;
		w(i * 4 + 3) = -c*VL / WL;

	}


	return w;
}
//////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////