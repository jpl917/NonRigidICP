#ifndef NONRIGID_ICP_H
#define NONRIGID_ICP_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <ctime>
#include <string>
#include <cmath>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SVD>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <Eigen/CholmodSupport>

#include <flann/flann.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/concept_check.hpp>

#include "TriMesh.h"
#include "TriMesh_algo.h"

using namespace std;
using namespace trimesh;

#define VERBOSE

#define NUM_LANDMARKS 51

#define THRESHOLD_DIST 20     //less
#define THRESHOLD_NORM 0.5     //greater

class MyTimer
{
public:
    timespec t0, t1; 
    MyTimer() {}
    double time_ms;
    double time_s;
    void start() {
		clock_gettime(CLOCK_REALTIME, &t0);
    }
    void end() {
		clock_gettime(CLOCK_REALTIME, &t1);
		time_ms = t1.tv_sec * 1000 + t1.tv_nsec/1000000.0 - (t0.tv_sec * 1000 + t0.tv_nsec/1000000.0);
		time_s = time_ms/1000;
		cout<<"   Time: "<<time_s<<" s"<<endl;
    }
};


MyTimer mytimer;

class NonrigidICP
{
private:
	trimesh::TriMesh* 					m_src;
	trimesh::TriMesh*              	 	m_dst;
	
	
	int 								m_vertsNum;						
	std::vector<trimesh::TriMesh::Face> m_faces;
	std::vector<std::pair<int,int> >    m_edges;
	
	
	//cv::Mat 							m_image;
	//Eigen::Matrix<double,3,4> 		m_projMatrix;
	
	//for landmark term
	std::vector<int>					Dl_Idx;  //template
	std::vector<Eigen::Vector3d> 		Dl;  
	
	std::vector<int>					Ul_Idx;
	std::vector<Eigen::Vector3d> 		Ul;  //target
	
	
	std::vector<float>					m_weights;
	std::vector<std::pair<int,int> >    m_soft_corres;
	
	
	double* 							m_dataset;
	flann::Index<flann::L2<double> >* 	m_kd_flann_index;
	
	//for rough ICP
	Eigen::Matrix3d 					ini_R;
	Eigen::Vector3d 					ini_t;
	
public:
	NonrigidICP(trimesh::TriMesh* _src, trimesh::TriMesh* _dst,
				std::vector<int> _Dl_Idx, std::vector<Eigen::Vector3d>  _Dl,
				std::vector<int> _Ul_Idx, std::vector<Eigen::Vector3d>  _Ul)
	{
		m_src = _src;
		m_dst = _dst;
		
		m_src->clear_colors();
		
		//m_image=_image;
		//m_projMatrix=_projMatrix;
		
		Dl_Idx = _Dl_Idx;
		Dl = _Dl;
		
		Ul_Idx = _Ul_Idx;
		Ul = _Ul;

		m_vertsNum = m_src->vertices.size();
		m_faces=m_src->faces;
		
		for(size_t i=0; i<m_faces.size();i++)
		{
			trimesh::TriMesh::Face face=m_faces[i];
			int a=face[0];
			int b=face[1];
			int c=face[2];
			
			m_edges.push_back(std::pair<int,int>(a,b));
			m_edges.push_back(std::pair<int,int>(b,c));
			m_edges.push_back(std::pair<int,int>(c,a));
		}
		
		buildKdTree(m_dst);
	}
	
	~NonrigidICP()
	{
		releaseKdTree();
	}
	
	
	void buildKdTree(trimesh::TriMesh* mesh)
	{
		std::vector<trimesh::point> verts = mesh->vertices;
		
		int vertNum=verts.size();
		m_dataset= new double[3*vertNum];
		
		for(int i=0; i<vertNum; i++)
		{
			m_dataset[3*i]=verts[i][0];
			m_dataset[3*i+1]=verts[i][1];
			m_dataset[3*i+2]=verts[i][2];
		}
		
		flann::Matrix<double> flann_dataset(m_dataset,vertNum,3);
		m_kd_flann_index=new flann::Index<flann::L2<double> >(flann_dataset,flann::KDTreeIndexParams(1));
		m_kd_flann_index->buildIndex();
		return;
	}
	
	void releaseKdTree()
	{
		delete[] m_dataset;
		delete m_kd_flann_index;
		m_kd_flann_index=NULL;
	}
	
		
	void resizeSrcMesh()
	{
		
		Eigen::Vector3d mean_Ul(0,0,0);
		Eigen::Vector3d mean_Dl(0,0,0);
		for(size_t i=0;i<Ul.size();i++)
		{
			mean_Ul += Ul[i];
			mean_Dl += Dl[i];
		}
		
		mean_Ul /= Ul.size();
		mean_Dl /= Dl.size();
		
		Eigen::Vector3d AA,BB;
		double meanRadius_AA=0, meanRadius_BB=0;
		for(size_t i=0; i<Ul.size();i++)
		{
			AA = Ul[i]-mean_Ul;
			BB = Dl[i]-mean_Dl;
			
			meanRadius_AA += sqrt(pow(AA(0),2)+pow(AA(1),2)+pow(AA(2),2));
			meanRadius_BB += sqrt(pow(BB(0),2)+pow(BB(1),2)+pow(BB(2),2));
		}
		
		float ratio = meanRadius_AA/meanRadius_BB; 
		cout<<"Scale ratio: "<<ratio<<endl;
		
		trimesh::scale(m_src, ratio);
		
		//Dl template
		for(size_t i=0; i<Dl_Idx.size(); i++)
		{
			int idx = Dl_Idx[i];
			trimesh::point xyz= m_src->vertices[idx];
			Dl[i]=Eigen::Vector3d(xyz[0],xyz[1],xyz[2]);
		}
	}
	
	

	void initialRigidICP()
	{
		Eigen::Vector3d mean_Ul(0,0,0);
		Eigen::Vector3d mean_Dl(0,0,0);
		for(size_t i=0;i<Ul.size();i++)
		{
			mean_Ul += Ul[i];
			mean_Dl += Dl[i];
		}
		
		mean_Ul /= Ul.size();
		mean_Dl /= Dl.size();
		
		//cout<<mean_Dl[0]<<" "<<mean_Dl[1]<<" "<<mean_Dl[2]<<endl;
		Eigen::Matrix<double,NUM_LANDMARKS,3> AA,BB;
		for(size_t i=0; i<Ul.size();i++)
		{
			AA(i,0) = Ul[i](0)-mean_Ul[0];
			AA(i,1) = Ul[i](1)-mean_Ul[1];
			AA(i,2) = Ul[i](2)-mean_Ul[2];
			
			BB(i,0) = Dl[i](0)-mean_Dl[0];
			BB(i,1) = Dl[i](1)-mean_Dl[1];
			BB(i,2) = Dl[i](2)-mean_Dl[2];
		}
		
		Eigen::Matrix3d W = Eigen::Matrix3d::Zero();  //= AA.transpose()*BB;
		for(size_t i=0; i<Ul.size(); i++)
		{
			W += Eigen::Vector3d(AA(i,0),AA(i,1),AA(i,2)) * Eigen::Vector3d(BB(i,0),BB(i,1),BB(i,2)).transpose();
		}
		
		Eigen::FullPivLU<Eigen::Matrix3d> lu_decomp(W);	
		Eigen::JacobiSVD<Eigen::Matrix3d> svd(W, Eigen::ComputeFullU| Eigen::ComputeFullV);
		
		Eigen::Matrix3d U = svd.matrixU();
		Eigen::Matrix3d V = svd.matrixV();
		
		ini_R = U * V.transpose();
		ini_t = -ini_R*mean_Dl+mean_Ul;
		
		if(ini_R.determinant()<0)
		{
			for(int j=0; j<3;j++) U(2,j) *= -1;
			ini_R = U*V.transpose();
		}
		
		cout<<"Initial Rigid ICP finished"<<endl;
	}

	
	void init()
	{
		//loadLandmarks();
		resizeSrcMesh();
		initialRigidICP();

		
#ifdef VERBOSE		
		ofstream fout("InitialAlign.obj");
#endif		
		
		
		for(int i=0;i<m_vertsNum;i++)
		{
			trimesh::point xyz= m_src->vertices[i];
			Eigen::Vector3d point(xyz[0], xyz[1], xyz[2]);
			Eigen::Vector3d icp_result = ini_R * point + ini_t;
			
			for(int d=0; d<3; d++)
				m_src->vertices[i][d] = icp_result(d);

#ifdef VERBOSE	
			fout<<"v "<<icp_result[0]<<" "<<icp_result[1]<<" "<<icp_result[2]<<endl;
#endif			
			
		}
		
#ifdef VERBOSE	
		fout.close();
#endif
		
		
		for(size_t i=0; i<Dl.size(); i++)
		{
			Eigen::Vector3d point(Dl[i](0),Dl[i](1),Dl[i](2));
			Eigen::Vector3d icp_result = ini_R * point + ini_t;
			
			Dl[i](0)=icp_result[0];
			Dl[i](1)=icp_result[1];
			Dl[i](2)=icp_result[2];
		}
	}

	void getKdCorrespondences()
	{
		m_weights.resize(m_vertsNum);
		for (size_t i = 0; i < m_weights.size(); ++i) m_weights[i] = 1.0f;
		
		const int knn=1;
		flann::Matrix<double> query(new double[3],1,3);
		flann::Matrix<int>    indices(new int[query.rows*knn],query.rows,knn);
		flann::Matrix<double> dists(new double[query.rows*knn],query.rows,knn);
	
		m_soft_corres.clear();
		m_soft_corres.resize(m_src->vertices.size());
		
		m_src->need_normals();
		m_dst->need_normals();
		
		for(size_t i=0; i<m_src->vertices.size();i++)
		{
			for(int j=0; j<3;j++)
				query[0][j]=m_src->vertices[i][j];
			
			m_kd_flann_index->knnSearch(query,indices,dists,1,flann::SearchParams(flann::FLANN_CHECKS_UNLIMITED));
			
			m_soft_corres[i]=std::make_pair(i,indices[0][0]);

			Eigen::Vector3d n1(m_src->normals[i][0],m_src->normals[i][1],m_src->normals[i][2]);
			Eigen::Vector3d n2(m_dst->normals[indices[0][0]][0],m_dst->normals[indices[0][0]][1],m_dst->normals[indices[0][0]][2]);
		
		
			//if( dists[0][0] > THRESHOLD_DIST  || abs(n1.dot(n2))<=THRESHOLD_NORM )  
			//	m_weights[i]=0.0f;
			
			
			//if(normal1.dot(normal2)<0) m_src->normals[i]*=-1;

			if( m_dst->is_bdy(indices[0][0]) || m_src->is_bdy(i))  m_weights[i]=0.0f;
			
		}
		
		delete[] query.ptr();
		delete[] indices.ptr();
		delete[] dists.ptr();
		
		double sum=0;
		for (size_t i = 0; i < m_weights.size(); ++i) sum+=m_weights[i];
			
		cout<<"Find matches:"<<sum<<endl;
		return;
	}


	// To do nonrigid icp registration
	Eigen::MatrixX3d compute(double alpha, double beta, double gamma)
	{
// 		m_faces.clear();
// 		m_edges.clear();
// 		m_faces=m_src->faces;
// 		
// 		for(size_t i=0; i<m_faces.size();i++)
// 		{
// 			trimesh::TriMesh::Face face=m_faces[i];
// 			int a=face[0];
// 			int b=face[1];
// 			int c=face[2];
// 			
// 			m_edges.push_back(std::pair<int,int>(a,b));
// 			m_edges.push_back(std::pair<int,int>(b,c));
// 			m_edges.push_back(std::pair<int,int>(c,a));
// 		}
		
		
		
		int n = m_src->vertices.size();
		int m = m_edges.size();
		int l = 0;
		if(Dl.size()>0 && Dl.size()==Ul.size())  l=Dl.size();
		
		
		Eigen::MatrixX3d X(4*n,3);
		X.setZero();
		
		
		bool loop = true;
		int  iter = 0;
		while(loop)
		{
			
			mytimer.start();
			getKdCorrespondences();
			mytimer.end();
			
			
			mytimer.start();
			cout<<"Construct AtA and Atb"<<endl;
			
			Eigen::SparseMatrix<double> A(4*m + n + l, 4*n);

			//1.alpha_M_G
			std::vector< Eigen::Triplet<double> > alpha_M_G;
			for (int i = 0; i < m; ++i)
			{
				int a = m_edges[i].first;
				int b = m_edges[i].second;

				for (int j = 0; j < 3; j++) 
				{
					alpha_M_G.push_back(Eigen::Triplet<double>(i*4 + j, a*4 + j, alpha));
					alpha_M_G.push_back(Eigen::Triplet<double>(i*4 + j, b*4 + j, -alpha));
				}
				
				alpha_M_G.push_back(Eigen::Triplet<double>(i*4 + 3, a*4 + 3, alpha * gamma));	
				alpha_M_G.push_back(Eigen::Triplet<double>(i*4 + 3, b*4 + 3, -alpha * gamma));
			}


			//2.W_D
			std::vector< Eigen::Triplet<double> > W_D;
			for (int i = 0; i < n; ++i)
			{				
				trimesh::point xyz = m_src->vertices[i];

				double weight = m_weights[i];
				if(weight==0)  weight=1;

				for (int j = 0; j < 3; ++j) 
					W_D.push_back(Eigen::Triplet<double>(4*m + i, i*4 + j, weight * xyz[j]));
				
				W_D.push_back(Eigen::Triplet<double>(4*m + i, i*4 + 3, weight*1));
			}

			
			//3.beta_D_L
			std::vector<Eigen::Triplet<double> > beta_D_L;
			for(int i = 0; i < l; i++)
			{
				for(int j = 0; j < 3; j++)
					beta_D_L.push_back(Eigen::Triplet<double>(4*m+n+i, Dl_Idx[i]*4+j, beta*Dl[i](j)));
				
				beta_D_L.push_back(Eigen::Triplet<double>(4*m+n+i, Dl_Idx[i]*4+3, beta));
			}
			

			std::vector< Eigen::Triplet<double> > _A = alpha_M_G;
			_A.insert(_A.end(), W_D.begin(), W_D.end());
			_A.insert(_A.end(), beta_D_L.begin(),beta_D_L.end());
			A.setFromTriplets(_A.begin(), _A.end());
			

			//for the B
			Eigen::MatrixX3d B = Eigen::MatrixX3d::Zero(4*m + n + l, 3);
			for (int i = 0; i < n; ++i)
			{
				
				int idx = 0;
				trimesh::point xyz;

				double weight = m_weights[i];
				if(weight != 0)
				{
					idx = m_soft_corres[i].second;
					xyz = m_dst->vertices[idx];
				}
				else
				{
					weight = 1;
					idx = m_soft_corres[i].first;
					xyz = m_src->vertices[idx];
				}
				
				for (int j = 0; j < 3; j++)  
					B(4*m + i, j) = weight * xyz[j];
			}
			
			for(int i = 0; i < l; i++)
			{
				for (int j = 0; j < 3; j++) 
					B(4*m + n + i, j) = beta * Ul[i](j);
			}
		
			Eigen::SparseMatrix<double> AtA = Eigen::SparseMatrix<double>(A.transpose()) * A;
			Eigen::MatrixX3d AtB = Eigen::SparseMatrix<double>(A.transpose()) * B;

		
			//Eigen::ConjugateGradient< Eigen::SparseMatrix<double> > solver;
			//Eigen::SimplicialCholesky<Eigen::SparseMatrix<double> > solver;
			Eigen::CholmodSupernodalLLT< Eigen::SparseMatrix<double> > solver;
			solver.compute(AtA);

			mytimer.end();
			
			
			//temporal X
			Eigen::MatrixX3d TmpX(4*n,3);
			TmpX=X;
			

			cout<<"Solve Ax=b  "<<endl;
			mytimer.start();
			X = solver.solve(AtB);
			mytimer.end();


			Eigen::Matrix3Xd Xt = X.transpose();
			for (int i = 0; i < n; ++i)
			{
				trimesh::point xyz=m_src->vertices[i];
				
				Eigen::Vector4d point(xyz[0], xyz[1], xyz[2], 1.0f);
				Eigen::Vector3d result = Xt.block<3, 4>(0, 4*i) * point;
			
				for(int d=0; d<3; d++)
				{
					m_src->vertices[i][d]=result[d];
				}
				
				
			}
			
			
			for(size_t i=0; i<Dl.size(); i++)
			{
				int idx = Dl_Idx[i];
				trimesh::point xyz= m_src->vertices[idx];
				Dl[i]=Eigen::Vector3d(xyz[0],xyz[1],xyz[2]);
			
			}
			
			cout<<"X Change:"<<(X-TmpX).norm()<<endl<<endl;
			if((X-TmpX).norm()<2)loop=0;
			
			iter++;
		}
		
		//cout<<"Sub iter:"<<iter<<endl;
		return X;
	}

};

#endif

