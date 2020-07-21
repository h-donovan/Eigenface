#include<iostream>
#include<dirent.h>
#include<sys/types.h>
#include<string>
#include<fstream>
#include<vector>
#include<math.h>
#include<list>

#include"jacobi.h"
#include"Matrix/Matrix.h"
#include"pgm/Image_PGM.h"


//SAMPLE SIZE is the number of images used to generate your model
//IMG_SIZE is the equal to the width*height of each image in your model.

//This requires .pgm filetype and some functions were designed specifically for an assignment and don't make sense outside of the assignment parameters.

#define SAMPLE_SIZE 1203
#define IMG_SIZE 320




using namespace std;

void list_dir(const char*);
void read_Face(const char*, Matrix&);//computes average face
void add_Face(const char*, Matrix*, int);//takes file location, Matrix array, and size
void MatToImg(ImageType&, Matrix&);
void xFace(const char*, Matrix*, int, string*);
void MtI_reconstruct(ImageType& img, Matrix& Mat, const int row, const int col);
void normalize(Matrix& vec);

//vector<int> get_dir(const char* path);

void start(Matrix* x, Matrix& x_bar, Matrix* Phi, int dims, string* id, const char* source);

void writeA(Matrix& A, Matrix* Phi, Matrix& ATA, const char* A_name, const char* ATA_name);
void readA(Matrix& A, Matrix& ATA, const char* A_name, const char* ATA_name);

void writeU(Matrix& A, Matrix &ATA, Matrix* u, double* lambda, int dims, const char* u_name, const char* lambda_name);
void readU(Matrix* u, double* lambda, int dims, int size,  const char* u_name, const char* lambda_name);

void writeOmega(Matrix* Omega, Matrix* Phi, Matrix *u, int dims, int size, const char* Omega_name);
void readOmega(Matrix* Omega, int dims, int size, const char* Omega_name);

double Euc_dist(Matrix& A, Matrix& B);

int topPercent(double percent, double* lambda, int size);

void project(Matrix& image, Matrix& Omega, Matrix* u, int k);

double Mah_dist(Matrix& A, Matrix& B, double* lam, int k);

struct val_and_index{
	double val;
	int index;
};

bool operator<(const val_and_index& a, const val_and_index& b){return a.val < b.val;}


int main(){

//Sample:
	cout<<"Start..."<<endl;
	string* x_id = new string[SAMPLE_SIZE];
	Matrix* x = new Matrix[SAMPLE_SIZE];
	Matrix* Phi = new Matrix[SAMPLE_SIZE];
	Matrix x_bar(IMG_SIZE,1);

	start(x, x_bar, Phi, SAMPLE_SIZE, x_id, "Faces_FA_FB/fa_L/");//initiates all images into x matrix, determines the mean, and generates x-x_bar to Phi.  Saves to file location.


	Matrix A(IMG_SIZE, SAMPLE_SIZE);
	Matrix AtA(IMG_SIZE, SAMPLE_SIZE);

	writeA(A, Phi, AtA, "d/A.bin", "d/ATA.bin");//will write to A and A.Transpose * A using Phi
//	readA(A, AtA, "d/A.bin", "d/ATA.bin");//will read from file location to Matrices.  writeA takes a long time so it's worth saving the data externally


	Matrix* u  = new Matrix[SAMPLE_SIZE];
	double* lambda = new double[SAMPLE_SIZE+1];//saves eigenvector coefficients.  Jacobi requires n+1 samples.
	writeU(A, AtA, u, lambda, SAMPLE_SIZE, "d/u.bin", "d/lambda.bin");//uses A and AtA to calculate u (mu) and the eigen coefficients.  Saves to file locations.
//	readU(u, lambda, 1203, 320, "d/u.bin", "d/lambda.bin");

	Matrix* Omega = new Matrix[SAMPLE_SIZE];
//	writeOmega(Omega, Phi, u, 1203, 1203, "d/Omega.bin");
	readOmega(Omega, SAMPLE_SIZE, SAMPLE_SIZE, "d/Omega.bin");

	cout<<u[0].arr[0][0]<<endl;



//The following tests are specific to the dataset and requirements used in this assignment.
//As such, they aren't super relevant to a lot of the data you may be looking for.

	for(int i=0; i<SAMPLE_SIZE; i++){
		double min=1000;
		double max=-1000;
		for(int j=0; j<320; j++){
			if(u[i].arr[j][0] > max){
				max = u[i].arr[j][0];
			}
			else if(u[i].arr[j][0] < min){
				min=u[i].arr[j][0];
			}
		}
		for(int j=0; j<320; j++){
			u[i].arr[j][0] = (u[i].arr[j][0]-min)/(max-min)*255;
		}
	}






	for(int i=0; i<10; i++){
		char num = i+'0';
		string best = "L_best";
		string worst = "L_worst";
		best += num;
		worst += num;
		best += ".pgm";
		worst += ".pgm";

		ImageType good(20,16,255);//the 20,16 is the image size of the sample i used.  that's why IMG_SIZE is 320.  This was a test on very small images.  you would change them to yours.
		ImageType bad(20,16,255);

		MtI_reconstruct(good, u[i], 20, 16);
		MtI_reconstruct(bad, u[1193+i], 20, 16);

		writeImage(best.c_str(), good);
		writeImage(worst.c_str(), bad);
		
	}


	delete[] x;
	delete[] x_id;
	delete[] Phi;
	delete[] u;
	delete[] lambda;
	delete[] Omega;


	cout<<"Ending main."<<endl;
	return 0;


}





double Mah_dist(Matrix& A, Matrix& B, double* lam, int k)//Calculates Mahalanobis distance between A and B
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


int topPercent(double percent, double* lambda, int size){
	double sum=0.0;
	for(int i=1; i<size+1; i++){
		sum += lambda[i];
	}
	double numerator=0.0;
	int k;
	for(k=1; k<size+1; k++){
		numerator += lambda[k];
		if(numerator/sum > percent) return k;
	}
	return -1;
}

void project(Matrix& image, Matrix& Omega, Matrix* u, int k){
	for(int i=0; i<k; i++){
		image = image + (Omega.arr[i][0] * u[i]);
	}
}


double Euc_dist(Matrix& A, Matrix& B){//calculates Euclidean distance between 2 matrices
	double sum=0.0;
	int row=A.row;
	int col = A.col;
	Matrix temp(row, col);//assumes they're the same size
	temp = A-B;
	for(int i=0; i<row; i++){
		for(int j=0; j<col; j++){
			sum+=temp.arr[i][j]*temp.arr[i][j];
		}
	}
	return sqrt(sum);
}


void writeOmega(Matrix* Omega, Matrix* Phi, Matrix *u, int dims, int size, const char* Omega_name){
	//dims = 1153
	//size = 1196
	cout<<"Calculating Omega..."<<endl;
	for(int i=0; i<size; i++){
		Omega[i] = Matrix(dims, 1);
	}
	for(int i=0; i<size; i++){
		for(int j=0; j<dims; j++){
			Omega[i].arr[j][0] = (u[j].T() * Phi[i]).arr[0][0];
		}
		cout<<i<<endl;
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
	int size = A.row;//2880
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

	cout<<"Jacobi dim="<<dims<<"..."<<endl;
	jacobi(ATA_, dims, lambda, V);

	cout<<"Calculating u..."<<endl;

	Matrix* vi = new Matrix[dims];
	for(int i=1; i<dims+1; i++){
		cout<<i<<endl;
		vi[i-1] = Matrix(dims, 1);
		for(int j=1; j<dims+1; j++){
			vi[i-1].arr[j-1][0] = V[j][i];
		}
	}

	for(int i=0; i<dims; i++){
		u[i] = Matrix(size, 1);
		u[i] = A*vi[i];//2880x1203 * 1203x1 = 2880x1
	}

	cout<<"Normalizing u..."<<endl;
	for(int i=0; i<dims; i++){
		normalize(u[i]);
	}


	ofstream u_out;
	ofstream lambda_out;

	u_out.open(u_name, ios::out|ios::binary);
	lambda_out.open(lambda_name, ios::out|ios::binary);

	cout<<"Writing u..."<<endl;
	//gonna try writing horizontally since this is causing so much trouble
	//Write u.  u[1203] (2880,1)
	for(int i=0; i<dims; i++){
		for(int j=0; j<size; j++){
			u_out << u[i].arr[j][0] << '\t';
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
	delete[] vi;
}

void readU(Matrix* u, double* lambda, int dims, int size, const char* u_name, const char* lambda_name){
	ifstream u_in;
	ifstream lambda_in;

	u_in.open(u_name, ios::in|ios::binary);
	lambda_in.open(lambda_name, ios::in|ios::binary);

	for(int i=0; i<dims+1; i++){
		lambda_in >> lambda[i];
	}

	//remember to read horizontally
	for(int i=0; i<dims; i++){
		u[i] = Matrix(size, 1);
	}
	for(int i=0; i<dims; i++){
		for(int j=0; j<size; j++){
			u_in >> u[i].arr[j][0];
		}
	}

//	ATA_in.close();
	u_in.close();
	lambda_in.close();
}











void start(Matrix* x, Matrix& x_bar, Matrix* Phi, int dims, string *id, const char* source){
	cout<<"Start..."<<endl;
	xFace(source, x, dims, id);
	
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
	int size=A.row;//2880		2880

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
	int size = A.row;//2880

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


void normalize(Matrix& vec){//assumes it's vector eg, vec = Matrix(N, 1).  Normalizes a matrix
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

void xFace(const char *path, Matrix* x, int size, string* id){
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
			id[ind] = location.substr(location.find_last_of('/')+1, 5);
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
