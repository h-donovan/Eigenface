#include<iostream>
#include "dirent.h"
#include<sys/types.h>
#include<string>
#include<fstream>
#include<vector>
#include<math.h>

#include"jacobi.h"
#include"Matrix.h"
#include"Image_PGM.h"
#include "ProgressBar.h"

#define FA_SIZE 1204
#define FA2_HSIZE 1138
#define FB_SIZE 1195

#define HR_ARRAY_SIZE 2880
#define LR_ARRAY_SIZE 320


#define FA_FOLDER "Faces_FA_FB/fa_H/"
#define FAL_FOLDER "Faces_FA_FB/fa_L/"
#define FB_FOLDER "Faces_FA_FB/fb_H/"
#define FBL_FOLDER "Faces_FA_FB/fb_L/"
#define FA2_FOLDER "Faces_FA_FB/fa2_H/"
#define FA2L_FOLDER "Faces_FA_FB/fa2_L/"

using namespace std;

void list_dir(const char*);
void read_Face(const char*, Matrix&);//computes average face
void add_Face(const char*, Matrix*, int);//takes file location, Matrix array, and size
void MatToImg(ImageType&, Matrix&);
void xFace(const char*, Matrix*, int);
void MtI_reconstruct(ImageType& img, Matrix& Mat, const int row, const int col);
void normalize(Matrix& vec);

void start(Matrix* x, Matrix& x_bar, Matrix* Phi, int dims, const char* source);

void writeA(Matrix& A, Matrix* Phi, Matrix& ATA, const char* A_name, const char* ATA_name);
void readA(Matrix& A, Matrix& ATA, const char* A_name, const char* ATA_name);

void writeU(Matrix& A, Matrix &ATA, Matrix* u, double* lambda, int dims, const char* u_name, const char* lambda_name);
void readU(Matrix* u, double* lambda, int dims, int size,  const char* u_name, const char* lambda_name);

void writeOmega(Matrix* Omega, Matrix* Phi, Matrix *u, int dims, int size, const char* Omega_name);
void readOmega(Matrix* Omega, int dims, int size, const char* Omega_name);

void Mah_dist_execution(Matrix* TestPhi, Matrix* HostPhi, Matrix& bar, Matrix* Omega, Matrix* u, double* lambda, int lambda_size,
	double low_threshold, double high_threshold, int arraySize, double eigen_percent, const char* folder, const char* output_file, 
	vector<int> intruder_list);
int calc_eigenCount(double percentage, double* lambda, int size);
bool isIntruder(int id, vector<int>knownIntruders);
vector<int> buildIntruderList(vector<int> original, vector<int> test_list);

double Euc_dist(Matrix& A, Matrix& B, int dims);
double Mah_dist_test(Matrix& image, Matrix& test, Matrix& mean, Matrix* Omega, Matrix* u, int k, double* lam, int dims);
double Euc_dist_test(Matrix& image, Matrix& test, Matrix& mean, Matrix* Omega, Matrix* u, int k, int dims);
double Mah_dist(Matrix& A, Matrix& B, double* lam, int dims);

void project(Matrix& image, Matrix& Omega, Matrix* u, int k);


vector<int> get_dir(const char* path);

int main(){


	Matrix *fb_H = new Matrix[FB_SIZE];//testing data
	Matrix fb_H_bar(HR_ARRAY_SIZE,1);
	Matrix* fb_H_Phi = new Matrix[FB_SIZE];
	start(fb_H, fb_H_bar, fb_H_Phi, FB_SIZE, FB_FOLDER);
	
	Matrix *fb_L = new Matrix[FB_SIZE];//testing data
	Matrix fb_L_bar(LR_ARRAY_SIZE,1);
	Matrix* fb_L_Phi = new Matrix[FB_SIZE];
	start(fb_L, fb_L_bar, fb_L_Phi, FB_SIZE, FBL_FOLDER);
	
/*	Matrix *fa_H = new Matrix[FA_SIZE];//testing data
	Matrix fa_H_bar(HR_ARRAY_SIZE,1);
	Matrix* fa_H_Phi = new Matrix[FA_SIZE];
	start(fa_H, fa_H_bar, fa_H_Phi, FA_SIZE, FB_FOLDER);
	
	Matrix *fa_L = new Matrix[FA_SIZE];//testing data
	Matrix fa_L_bar(LR_ARRAY_SIZE,1);
	Matrix* fa_L_Phi = new Matrix[FA_SIZE];
	start(fa_L, fa_L_bar, fa_L_Phi, FA_SIZE, FBL_FOLDER);


/*	Matrix *x = new Matrix[1203];
	Matrix x_bar(HR_ARRAY_SIZE,1);
	Matrix* Phi = new Matrix[1203];
	start(x, x_bar, Phi, 1203, "Faces_FA_FB/fa_H/");

	Matrix AtA(1203, 1203);
	Matrix A(HR_ARRAY_SIZE, 1203);
//	writeA(A, Phi, AtA, "a/A.bin", "a/ATA.bin");
	readA(A, AtA, "a/A.bin", "a/ATA.bin");

	Matrix* u = new Matrix[1203];
	double* lambda = new double[1204];
//	writeU(A, AtA, u, lambda, 1203,  "a/u.bin", "a/lambda.bin");
	readU(u, lambda, 1203, HR_ARRAY_SIZE, "a/u.bin", "a/lambda.bin");

	Matrix* Omega = new Matrix[1203];
//	writeOmega(Omega, Phi, u, 1203, 1203, "a/Omega.bin");	
	readOmega(Omega, 1203, 1203, "a/Omega.bin");*/

//====================================================================
//b at High Resolution

	cout<<endl<<"Part b:"<<endl;
	Matrix* fa2_H = new Matrix[FA2_HSIZE];
	Matrix* fa2_H_Phi = new Matrix[FA2_HSIZE];
	Matrix fa2_H_bar(HR_ARRAY_SIZE,1);

	start(fa2_H, fa2_H_bar, fa2_H_Phi, FA2_HSIZE, FA2_FOLDER);

	Matrix A_b(HR_ARRAY_SIZE, FA2_HSIZE);
	Matrix AtA_b(FA2_HSIZE, FA2_HSIZE);

//	writeA(A_b, fa2_H_Phi, AtA_b, "b/A.bin", "b/ATA.bin");
	readA(A_b, AtA_b, "b/A.bin", "b/ATA.bin");

	Matrix* u_b = new Matrix[FA2_HSIZE];
	double* lambda_b = new double[FA2_HSIZE+1];
//	writeU(A_b, AtA_b, u_b, lambda_b, FA2_HSIZE, "b/u.bin", "b/lambda.bin");
	readU(u_b, lambda_b, FA2_HSIZE, HR_ARRAY_SIZE, "b/u.bin", "b/lambda.bin");

	Matrix* Omega_b = new Matrix[FA2_HSIZE];
//	writeOmega(Omega_b, fa2_H_Phi, u_b, FA2_HSIZE, FA2_HSIZE, "b/Omega.bin");
//	writeOmega(Omega_testing_b, fb_H_Phi, u_b, 1153, FB_SIZE, "b/Omega_testing.bin");//int1 dims, int2 size,...,		probably wrong...

	readOmega(Omega_b, FA2_HSIZE, FA2_HSIZE, "b/Omega.bin");

// ==========================================================================
// b at Low Resolution
	
	cout << endl << "Part b low-res:" << endl;
	Matrix* fa2_L = new Matrix[FA2_HSIZE];
	Matrix* fa2_L_Phi = new Matrix[FA2_HSIZE];
	Matrix fa2_L_bar(LR_ARRAY_SIZE,1);

	start(fa2_L, fa2_L_bar, fa2_L_Phi, FA2_HSIZE, FA2L_FOLDER);

	Matrix A_bl(LR_ARRAY_SIZE, FA2_HSIZE);
	Matrix AtA_bl(FA2_HSIZE, FA2_HSIZE);

//	writeA(A_bl, fa2_L_Phi, AtA_bl, "b/A_l.bin", "b/ATA_l.bin");
	readA(A_bl, AtA_bl, "b/A_l.bin", "b/ATA_l.bin");

	Matrix* u_bl = new Matrix[FA2_HSIZE];
	double* lambda_bl = new double[FA2_HSIZE+1];
//	writeU(A_bl, AtA_bl, u_bl, lambda_bl, FA2_HSIZE, "b/u_l.bin", "b/lambda_l.bin");
	readU(u_bl, lambda_bl, FA2_HSIZE, LR_ARRAY_SIZE, "b/u_l.bin", "b/lambda_l.bin");

	Matrix* Omega_bl = new Matrix[FA2_HSIZE];
//	writeOmega(Omega_bl, fa2_L_Phi, u_bl, FA2_HSIZE, FA2_HSIZE, "b/Omega_l.bin");
//	writeOmega(Omega_testing_b, fb_H_Phi, u_b, 1153, FB_SIZE, "b/Omega_testing_l.bin");//int1 dims, int2 size,...,		probably wrong...

	readOmega(Omega_bl, FA2_HSIZE, FA2_HSIZE, "b/Omega_l.bin");

//113 elements needed for 95%
//do go talk to Bebis about this, but ed vs er???
//ek=||Phi - Phi_hat||
//if ek<T, then this person IS on the list
//Max Threshold = 7961.03
//Take an image,

	Matrix* Testing_bl_Phi = new Matrix[FB_SIZE];
	Matrix* Testing_b_Phi = new Matrix[FB_SIZE];
	for(int i=0; i<FB_SIZE; i++){
		Testing_b_Phi[i] = Matrix(HR_ARRAY_SIZE, 1);
		Testing_b_Phi[i] = fb_H[i] - fa2_H_bar;
		Testing_bl_Phi[i] = Matrix(LR_ARRAY_SIZE, 1);
		Testing_bl_Phi[i] = fb_L[i] - fa2_L_bar;
	}

	Matrix* Testing_b_Omega = new Matrix[FB_SIZE];
//	writeOmega(Testing_b_Omega, Testing_b_Phi, u_b, FA2_HSIZE, FB_SIZE, "b/Omega_testing.bin");
	readOmega(Testing_b_Omega, FA2_HSIZE, FB_SIZE, "b/Omega_testing.bin");
	Matrix* Testing_bl_Omega = new Matrix[FB_SIZE];
//	writeOmega(Testing_bl_Omega, Testing_bl_Phi, u_bl, FA2_HSIZE, FB_SIZE, "b/Omega_testing_l.bin");
	readOmega(Testing_bl_Omega, FA2_HSIZE, FB_SIZE, "b/Omega_testing_l.bin");

	cout << "Begining Mahala test." << endl;

	// threshold range: 0.0296 to 1.2058, approx.
	vector<int> a_ids = get_dir(FA_FOLDER);
	vector<int> b_ids = get_dir(FB_FOLDER);
	vector<int> a2_ids = get_dir(FA2_FOLDER);


	cout << "Building intruder lists." << endl;
	vector<int> intruders_ab = buildIntruderList(a_ids, b_ids);
	vector<int> intruders_a2b = buildIntruderList(a2_ids, b_ids);

	for (int i = 0; i <= 50; i++)
	{
		intruders_a2b.push_back(i);
	}

	Mah_dist_execution(Testing_b_Phi, fa2_H_Phi, fa2_H_bar, Testing_b_Omega, u_b, lambda_b, FA2_HSIZE + 1,  0.0, 0.55, HR_ARRAY_SIZE, 0.95, FB_FOLDER, "Maha H Test", intruders_a2b);
	Mah_dist_execution(Testing_bl_Phi, fa2_L_Phi, fa2_L_bar, Testing_bl_Omega, u_bl, lambda_bl, FA2_HSIZE + 1, 0.08, 0.481, LR_ARRAY_SIZE, 0.95, FBL_FOLDER, "Maha L Test", intruders_a2b);

//	int k=113;

	/*cout<<"T\tTP\tFP"<<endl;
	for(int i=0; i<7961; i+=300){
		int TP=0;
		int FP=0;
		//Test for FP (True positives)
		for(int j=0; j<FB_SIZE; j++){
			Matrix temp(2800, 1);
			for(int k=0; k<113; k++){
				temp = temp + (Testing_b_Omega[j].arr[k][0] * u[k]);
			}
			if(Euc_dist(temp, Testing_b_Phi[j], 2800) < i){
				FP++;
			}
		}
		//Test for TP (False positives)
		for(int j=0; j<1153; j++){
			Matrix temp(2800, 1);
			for(int k=0; k<113; k++){
				temp = temp + (Omega_b[j].arr[k][0] * u[k]);
			}
			if(Euc_dist(temp, fa2_H_Phi[j], 2800) < i){
				TP++;
			}
		}
		cout<<i<<'\t'<<TP<<'\t'<<FP<<endl;
		if(FP>=1195 && TP >=1152) break;
	}*/




//Clean Up
	delete[] fb_H;
	delete[] fb_H_Phi;
	delete[] fb_L;
	delete[] fb_L_Phi;

	/*delete[] x;
	delete[] Phi;
	delete[] u;
	delete[] lambda;
	delete[] Omega;*/

	delete[] fa2_H;
	delete[] fa2_H_Phi;
	delete[] u_b;
	delete[] lambda_b;
	delete[] Omega_b;

	delete[] Testing_b_Phi;
	delete[] Testing_b_Omega;
	
	delete[] fa2_L;
	delete[] fa2_L_Phi;
	delete[] u_bl;
	delete[] lambda_bl;
	delete[] Omega_bl;

	delete[] Testing_bl_Phi;
	delete[] Testing_bl_Omega;

	cout<<"Ending main."<<endl;
	return 0;


}




void Mah_dist_execution(Matrix* TestPhi, Matrix* HostPhi, Matrix& bar, Matrix* Omega, Matrix* u, double* lambda, int lambda_size,
	double low_threshold, double high_threshold, int arraySize, double eigen_percent, const char* folder, const char* output_file, 
	vector<int> intruder_list)
{
	int eigenCount = calc_eigenCount(eigen_percent, lambda, lambda_size);

	ofstream fout;
	string filename = "b/";
	filename.append(output_file);
	filename.append("_");
	filename.append(std::to_string((int)(eigen_percent * 100)));
	filename += ".txt";
	fout.open(filename);

	vector<int> b_ids = get_dir(folder);

	

	int positives = 0;
	int negatives = 0;
	int fp = 0;
	int fn = 0;
	for (double t = low_threshold; t <= high_threshold; t += 0.0025)
	{
		positives = 0;
		negatives = 0;
		fp = 0;
		fn = 0;
		for (int j = 0; j < FB_SIZE; j++)
		{
			bool p = false;
			for (int i = 0; i < FA2_HSIZE; i++)
			{
				double test = Mah_dist_test(TestPhi[j], HostPhi[i], bar, Omega, u, eigenCount, lambda, HR_ARRAY_SIZE);
				if (test < t)
				{
					positives++;
					bool fa_ = isIntruder(b_ids[j], intruder_list);
					if (fa_)
						fp++;

					p = true;
					break;
				}
			}

			if (!p)
			{
				negatives++;
				bool fa_ = isIntruder(b_ids[j], intruder_list);
				if (fa_)
					fn++;
			}


			printDetailedProgressBar(j + 1, FB_SIZE);
		}

		fout << t << "\t" << positives << "\t" << fp << "\t" << negatives << "\t" << fn << endl;
	}

	fout.close();
}

int calc_eigenCount(double percentage, double* lambda, int size)
{
	double denom = 0;
	double numer = 0;

	for (int i = 1; i <= size; i++)
	{
		denom += lambda[i];
	}

	double perc = 0;
	int c = 0;
	while (c <= size && perc < percentage)
	{
		numer += lambda[++c];
		perc = numer / denom;
	}

	return c;
}

vector<int> buildIntruderList(vector<int> original, vector<int> test_list)
{
	vector<int> list;
	for (int i = 0; i < test_list.size(); i++)
	{
		for (int k = 0; k < original.size(); k++)
		{
			if (test_list[i] == original[k])
			{
				list.push_back(test_list[i]);
				break;
			}
		}
	}

	return list;
}

bool isIntruder(int id, vector<int>knownIntruders)
{
	if (knownIntruders.empty())
		return false;

	for (int k = 0; k < knownIntruders.size(); k++)
	{
		if (id == k)
			return true;
	}

	return false;
}


void project(Matrix& image, Matrix& Omega, Matrix* u, int k){
	for(int i=0; i<k; i++){
		image = image + (Omega.arr[i][0] * u[i]);
	}
}


double Euc_dist(Matrix& A, Matrix& B, int dims){
	double sum=0.0;
	double curr;
	for(int i=0; i<dims; i++){
		curr = A.arr[i][0]-B.arr[i][0];
		sum += curr*curr;
	}
	return sqrt(sum);
}

double Mah_dist_test(Matrix& image, Matrix& test, Matrix& mean, Matrix* Omega, Matrix* u, int k, double* lam, int dims)
{
	//project(image, *Omega, u, k);
	return Mah_dist(image, test, lam, k);
}

double Euc_dist_test(Matrix& image, Matrix& test, Matrix& mean, Matrix* Omega, Matrix* u, int k, int dims)
{
	//project(image, *Omega, u, k);
	return Euc_dist(image, test, k);
}

double Mah_dist(Matrix& A, Matrix& B, double* lam, int k)
{
	double sum = 0.0;
	double curr;
	for (int i = 0; i < k; i++)
	{
		curr = A.arr[i][0] - B.arr[i][0];
		sum += ((curr * curr) / lam[i+1]);
	}
	return sqrt(sum);
}






void writeOmega(Matrix* Omega, Matrix* Phi, Matrix *u, int dims, int size, const char* Omega_name){
	//dims = 1153
	//size = FB_SIZE
	cout<<"Calculating Omega..."<<endl;
	for(int i=0; i<size; i++){
		Omega[i] = Matrix(dims, 1);
	}
	for(int i=0; i<size; i++){
		for(int j=0; j<dims; j++){
			Omega[i].arr[j][0] = (u[j].T() * Phi[i]).arr[0][0];
		}
		printDetailedProgressBar(i + 1, size);
	}

	cout<<"Writing Omega..."<<endl;
	ofstream Omega_out;
	Omega_out.open(Omega_name, ios::out|ios::binary);

	for(int i=0; i<dims; i++){
		for(int j=0; j<size; j++){
			Omega_out << Omega[j].arr[i][0]<<'\t';
		}
		Omega_out<<endl;
	}
	Omega_out.close();
}

void readOmega(Matrix* Omega, int dims, int size, const char* Omega_name){
	cout<<"Reading Omega..."<<endl;
	ifstream Omega_in;
	Omega_in.open(Omega_name, ios::in|ios::binary);

	for(int i=0; i<size; i++){
		Omega[i] = Matrix(dims, 1);
	}
	for(int i=0; i<dims; i++){
		for(int j=0; j<size; j++){
			Omega_in >> Omega[j].arr[i][0];
		}
	}
	Omega_in.close();
}










void writeU(Matrix& A, Matrix &ATA, Matrix* u, double* lambda, int dims, const char* u_name, const char* lambda_name){
	int size = A.row;
	double** ATA_ = new double*[dims+1];
	double** V = new double*[dims+1];
	for(int i=0; i<dims+1; i++){
		ATA_[i] = new double[dims+1];
		V[i] = new double[dims+1];
	}

	for(int i=1; i<dims+1; i++){
		for(int j=1; j< dims+1; j++){
			ATA_[i][j] = ATA.arr[i-1][j-1];
		}
	}

	cout<<"Jacobi..."<<endl;
	jacobi(ATA_, dims, lambda, V);

	cout<<"Calculating u..."<<endl;
	Matrix v_i(dims, 1);
	for(int i=0; i<dims; i++){
		for(int j=0; j<dims; j++){
			v_i.arr[j][0] = V[j+1][i+1];
		}
		u[i] = Matrix(size, 1);
		u[i] = A*v_i;
	}

	cout<<"Normalizing u..."<<endl;
	for(int i=0; i<dims; i++){
		normalize(u[i]);
	}

//	ofstream ATA_out;
	ofstream u_out;
	ofstream lambda_out;

//	ATA_out.open(ATA_name, ios::out|ios::binary);
	u_out.open(u_name, ios::out|ios::binary);
	lambda_out.open(lambda_name, ios::out|ios::binary);

//	cout<<"Writing ATA..."<<endl;
//	//Write ATA
//	for(int i=0; i<dims; i++){
//		for(int j=0; j<dims; j++){
//			ATA_out << ATA.arr[i][j] << '\t';
//		}
//		ATA_out << endl;
//	}
	cout<<"Writing u..."<<endl;
	//Write u
	for(int i=0; i<A.row; i++){
		for(int j=0; j<dims; j++){
			u_out << u[j].arr[i][0] << '\t';
		}
		u_out << endl;
	}
	cout<<"Writing lambda..."<<endl;
	for(int i=0; i<dims+1; i++){
		lambda_out << lambda[i] << endl;
	}

//	ATA_out.close();
	u_out.close();
	lambda_out.close();

	cout<<"Cleaning..."<<endl;
	for(int i=0; i<dims+1; i++){
		delete[] ATA_[i];
		delete[] V[i];
	}
	delete[] ATA_;
	delete[] V;

}

void readU(Matrix* u, double* lambda, int dims, int size, const char* u_name, const char* lambda_name){
//	ifstream ATA_in;
	ifstream u_in;
	ifstream lambda_in;

//	ATA_in.open(ATA_name, ios::in|ios::binary);
	u_in.open(u_name, ios::in|ios::binary);
	lambda_in.open(lambda_name, ios::in|ios::binary);

//	for(int i=0; i<dims; i++){
//		for(int j=0; j<dims; j++){
//			ATA_in >> ATA.arr[i][j];
//		}
//	}
	for(int i=0; i<dims+1; i++){
		lambda_in >> lambda[i];
	}

	for(int i=0; i<dims; i++){
		u[i] = Matrix(size, 1);
	}
	for(int i=0; i<size; i++){
		for(int j=0; j<dims; j++){
			u_in >> u[j].arr[i][0];
		}
	}

//	ATA_in.close();
	u_in.close();
	lambda_in.close();
}











void start(Matrix* x, Matrix& x_bar, Matrix* Phi, int dims, const char* source){
	cout<<"Start..."<<endl;
	xFace(source, x, dims);
	
	int size = x_bar.row;
	cout<<"X bar..."<<endl;
	for(int i=0; i<dims; i++){
		x_bar = x_bar + x[i];
	}
	x_bar = (1.0/dims) * x_bar;

	cout<<"Phi..."<<endl;
	for(int i=0; i<dims; i++){
		Phi[i] = Matrix(size, 1);
		Phi[i] = x[i] - x_bar;
	}

}












void writeA(Matrix& A, Matrix* Phi, Matrix& ATA, const char* A_name, const char* ATA_name){
	int dims=A.col;//1203		1153
	int size=A.row;//HR_ARRAY_SIZE		HR_ARRAY_SIZE

	cout<<"A..."<<endl;
	for(int i=0; i<size; i++){
		for(int j=0; j<dims; j++){
			A.arr[i][j] = Phi[j].arr[i][0];
		}
	}

	cout<<"ATA..."<<endl;
	ATA = A.T() * A;

	ofstream A_out;
	ofstream ATA_out;

	A_out.open(A_name, ios::out|ios::binary);
	ATA_out.open(ATA_name, ios::out|ios::binary);

	cout<<"Writing A..."<<endl;
	for(int i=0; i<size; i++){
		for(int j=0; j<dims; j++){
			A_out << A.arr[i][j]<<'\t';
		}
		A_out << endl;
	}

	cout<<"Writing ATA..."<<endl;
	for(int i=0; i<dims; i++){
		for(int j=0; j<dims; j++){
			ATA_out << ATA.arr[i][j] << '\t';
		}
		ATA_out << endl;
	}

	A_out.close();
	ATA_out.close();
	cout<<"Finished writing A & ATA"<<endl;
}

void readA(Matrix& A, Matrix& ATA, const char* A_name, const char* ATA_name){

	int dims = A.col;//1203
	int size = A.row;//HR_ARRAY_SIZE

	ifstream A_in;
	ifstream ATA_in;

	A_in.open(A_name, ios::in|ios::binary);
	ATA_in.open(ATA_name, ios::in|ios::binary);

	for(int i=0; i<size; i++){
		for(int j=0; j<dims; j++){
			A_in >> A.arr[i][j];
		}
	}
	for(int i=0; i<dims; i++){
		for(int j=0; j<dims; j++){
			ATA_in >> ATA.arr[i][j];
		}
	}

	A_in.close();
	ATA_in.close();
}




















void normalize(Matrix& vec){//assumes it's vector eg, vec = Matrix(N, 1)
	double sum=0.0;
	double curr;
	for(int i=0; i<vec.row; i++){
		curr = vec.arr[i][0];
		sum += curr*curr;
	}
	sum = sqrt(sum);
	vec = (1.0/sum) * vec;
}



void MatToImg(ImageType& img, Matrix& Mat){
	int row,col,d;
	img.getImageInfo(row,col,d);
	if(row!=Mat.row || col!=Mat.col){
		img.setImageInfo(Mat.row, Mat.col, d);
	}
	for(int i=0; i<Mat.row; i++){
		for(int j=0; j<Mat.col; j++){
			img.setVal(i,j,Mat.arr[i][j]);
		}
	}
}

void MtI_reconstruct(ImageType& img, Matrix& Mat, const int row, const int col){
	img.setImageInfo(row, col, 255);
	for(int i=0; i<row*col; i++){
		img.setVal(i/col, i%col, Mat.arr[i][0]);
	}
}

void xFace(const char *path, Matrix* x, int size){
	struct dirent *entry;
	DIR *dir = opendir(path);
	if(dir==NULL){
		cout<<"File not found."<<endl;
		return;
	}
	
	string location;
	ImageType temp;
	int ind=0;

	while((entry = readdir(dir)) != NULL && ind<size){
		location=path;
		location += entry->d_name;
		if(readImage(location.c_str(), temp)){
			int row, col, d;
			temp.getImageInfo(row, col, d);
			x[ind] = Matrix(row*col,1);
			for(int i=0; i<row*col; i++){
				temp.getVal(i/col,i%col,d);
				x[ind].arr[i][0] = d;
			}
			ind++;
		}
	}
}

void add_Face(const char *path, Matrix *Phi, int size){
	struct dirent *entry;
	DIR *dir = opendir(path);
	if(dir==NULL){
		cout<<"File not found."<<endl;
		return;
	}
	
	string location;
	ImageType temp;
	int x=0;

	while((entry = readdir(dir)) != NULL && x<size){
		location=path;
		location += entry->d_name;
		if(readImage(location.c_str(), temp)){
			int row, col, d;
			temp.getImageInfo(row, col, d);
			for(int i=0; i<row; i++){
				for(int j=0; j<col; j++){
					temp.getVal(i,j,d);
					Phi[x].arr[i][j] += d;
				}
			}
			x++;
		}
	}
}

void read_Face(const char *path, Matrix &val){
	struct dirent *entry;
	DIR *dir = opendir(path);
	if(dir==NULL){
		cout<<"File not found."<<endl;
		return;
	}

	int size=0;
	ImageType temp;
	string location;
	while((entry = readdir(dir)) != NULL){
		location=path;
		location += entry->d_name;

		if(readImage(location.c_str(), temp)){
			int row, col, d;
			temp.getImageInfo(row, col, d);
			for(int i=0; i<row; i++){
				for(int j=0; j<col; j++){
					temp.getVal(i,j,d);
					val.arr[i][j] += d;
				}
			}
			size++;
		}
	}
	cout<<size<<endl;
	val = (1.0/size)*val;
}

void list_dir(const char *path){
	struct dirent *entry;
	DIR *dir = opendir(path);

	if(dir==NULL){
		return;
	}

	while((entry = readdir(dir)) != NULL){
		cout<<entry->d_name<<endl;
	}
	closedir(dir);
}

//average face = 1/M * Sigma_n_to_M:(face n)



vector<int> get_dir(const char* path) {
	struct dirent* entry;
	DIR* dir = opendir(path);

	vector<int> list;

	if (dir == NULL)
	{
		cout << "No directory found: " << path << endl;
		exit(1);
	}

	while ((entry = readdir(dir)) != NULL)
	{
		string temp(entry->d_name);
		size_t p = temp.find_last_of("/\\");
		string t_sub = temp.substr(p + 1);

		if (t_sub._Equal(".") || t_sub._Equal(".."))
			continue;

		size_t v = t_sub.find_first_of("_");

		list.push_back(std::stoi(t_sub.substr(0, v)));
	}
	closedir(dir);

	return list;
}