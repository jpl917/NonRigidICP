#include "nonRigidICP.h"
#include "landmarkDetector.h"


using namespace std;

Eigen::Matrix<double,3,4> loadProjMatrix(string filename)
{
	Eigen::Matrix<double,3,4> P;
	
	ifstream fin(filename);
	string tmp;
	for(int i=0;i<4; i++)
	{
		if(i==0) fin>>tmp;
		else
		{
			fin>>P(i-1,0)>>P(i-1,1)>>P(i-1,2)>>P(i-1,3);
		}
	}
	fin.close();
	return P;
}


std::vector<int> loadLandmarksIdx(string filename)
{
	
	std::vector<int> landmarksIdx;
	int idx;
	ifstream fin(filename);
	for(int i=0; i<51; i++)  //68-17
	{
		fin>>idx;
		landmarksIdx.push_back(idx);
	}
	fin.close();
	
	return landmarksIdx;
}
	

int main(int argc, char**argv)
{
	
	std::string rootPath = "../demo";
	std::string imageFilename = rootPath+"/00000004.jpg";
	std::string matrixPFilename = rootPath+"/00000004.txt";
	std::string savePath = "../result/";
	
	std::string srcFilename = "../template/template.ply";   
	std::string dstFilename = rootPath + "/poisson_mesh.ply";   
	std::string srcLandmarksFilename = "../template/source.txt";

	string command;
    	command = "mkdir -p ../result/";
	system(command.c_str());	

	MyTimer timer;
	timer.start();
	
	
	trimesh::TriMesh *src, *dst;
	src = trimesh::TriMesh::read(srcFilename);
	dst = trimesh::TriMesh::read(dstFilename);
	
	src->need_neighbors();
	src->need_normals();
	
	dst->need_neighbors();
	dst->need_normals();
	
	
	
	//for the DL

	
	std::vector<int> Dl_Idx = loadLandmarksIdx(srcLandmarksFilename);
	
	std::vector<Eigen::Vector3d> Dl;
	//Dl template
	for(size_t i=0; i<Dl_Idx.size(); i++)
	{
		int idx = Dl_Idx[i];
		trimesh::point xyz= src->vertices[idx];
		Dl.push_back(Eigen::Vector3d(xyz[0],xyz[1],xyz[2]));
	}
	
	//for the Ul
	cv::Mat image = cv::imread(imageFilename);
	Eigen::Matrix<double,3,4> projMatrix = loadProjMatrix(matrixPFilename);
	LandmarkDetector landmarkDetector(dst,image,projMatrix);
	std::vector<Eigen::Vector3d> Ul = landmarkDetector.getLandmarks3D();
	std::vector<int> Ul_Idx(Ul.size(),0);
	
	cout<<"Dl size:"<<Dl.size()<<endl;
	cout<<"Ul size:"<<Ul.size()<<endl;
	
	NonrigidICP nricp(src, dst, Dl_Idx, Dl, Ul_Idx, Ul);
	nricp.init();
	
	//return 0;

	double max_alpha=1000;
	double min_alpha=50;
	double beta = 10.0;
	double gamma = 1.0;
	int step = 10;

	Eigen::MatrixX3d X(4* src->vertices.size(),3);
	X.setZero();
	
	//trimesh::subdiv(src,SUBDIV_LOOP);
	
	for (int i = 1; i <= step; ++i)
	{
		double alpha = max_alpha-i*(max_alpha-min_alpha)/step;
		
		std::cout << "*********************************"<<endl;
		std::cout<<"Iteration:" <<i<<"  alpha:"<<alpha<<endl;
		
		//Eigen::MatrixX3d TmpX=X;
		
		X = nricp.compute(alpha, beta, gamma);
		
		//if((X-TmpX).norm()<epsilon)break;
		
		std::string filename=savePath+std::to_string(i)+".ply";
		src->write(filename);
	}
	
	
	trimesh::subdiv(src,SUBDIV_LOOP);
	trimesh::lmsmooth(src,5);
	src->clear_normals();
	src->clear_colors();
	std::string filename=savePath+"result.ply";
	src->write(filename);
	
	
	cout<<"All Cost time:";
	timer.end();
	
	//for the whole head model
// 	{
// 		std::string rootPath = "../demo3";
// 		std::string imageFilename = rootPath+"/2.jpg";
// 		std::string matrixPFilename = rootPath+"/00000002.txt";
// 		std::string savePath = rootPath+"/result44/";
// 		
// 		std::string dstFilename = rootPath + "/head-right.obj";    
// 		
// 		trimesh::TriMesh *dst = trimesh::TriMesh::read(dstFilename);
// 
// 		dst->need_neighbors();
// 		dst->need_normals();
// 		
// 		cv::Mat image = cv::imread(imageFilename);
// 		Eigen::Matrix<double,3,4> projMatrix = loadProjMatrix(matrixPFilename);
// 		LandmarkDetector landmarkDetector(dst,image,projMatrix);
// 		std::vector<Eigen::Vector3d> Dl = landmarkDetector.getLandmarks3D();
// 		std::vector<int> Dl_Idx(Dl.size(),0);
// 		
// 		NonrigidICP nricp(dst, src, Dl_Idx, Dl, Ul_Idx, Ul);
// 		nricp.init();
// 		
// 		//return 0;
// 
// 		double max_alpha=1000;
// 		double min_alpha=10;
// 		double beta = 0.0;
// 		double gamma = 1.0;
// 		//double epsilon=1e-4;
// 		int step = 10;
// 
// 		Eigen::MatrixX3d X(4* dst->vertices.size(),3);
// 		X.setZero();
// 		for (int i = 1; i <= step; ++i)
// 		{
// 			double alpha = max_alpha-i*(max_alpha-min_alpha)/step;
// 			
// 			std::cout << "*********************************"<<endl;
// 			std::cout<<"Iteration:" <<i<<"  alpha:"<<alpha<<endl;
// 			X = nricp.compute(alpha, beta, gamma);
// 			std::string filename=savePath+std::to_string(i)+".ply";
// 			dst->write(filename);
// 		}
// 		
// 	}
	
	return 0;
}


