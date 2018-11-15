#include <fstream>
#include <algorithm>
#include <cmath>
#include <boost/graph/graph_concepts.hpp>

#include <dlib/image_processing/frontal_face_detector.h>
#include <dlib/image_processing/render_face_detections.h>
#include <dlib/image_processing.h>
#include <dlib/gui_widgets.h>
#include <dlib/image_io.h>
#include <dlib/opencv.h>

#include "TriMesh.h"
#include "TriMesh_algo.h"

#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

using namespace std;
using namespace dlib;

#define VERBOSE
#define NUM_POINTS 68

//2D distance 
bool isNear(Eigen::Vector3d& landmark, Eigen::Vector3d& point, float threshold)
{
	double d= sqrt(pow(landmark(0)-point(0),2)+ pow(landmark(1)-point(1),2));
	if(d>threshold) return false;
	return true;
}

bool cmp(const std::pair< trimesh::point, int>& a, const std::pair<trimesh::point, int>& b)
{
	if(a.first[0]+a.first[1]+a.first[2] > b.first[0]+b.first[1]+b.first[2]) return true;
	else return false;
}



class LandmarkDetector
{
private:
	trimesh::TriMesh*					m_mesh;			
	cv::Mat 							m_image;
	Eigen::Matrix<double,3,4> 			m_projMatrix;   //K[R|t]
	
	std::vector<Eigen::Vector3d> 		m_landmarks2D;  //68 2D image landmarks
	std::vector<Eigen::Vector3d> 		m_landmarks3D;
	
	//std::vector<int>					m_landmarks3D_Idx;
	
public:
	LandmarkDetector(trimesh::TriMesh* mesh,
				cv::Mat image,  Eigen::Matrix<double,3,4> _projMatrix)
	{
		m_mesh = mesh;
		m_image= image;
		m_projMatrix = _projMatrix;
		
		detect2DLandmarks();
		detect3DLandmarks();
		
	}
	
	~LandmarkDetector(){}
	
	std::vector<Eigen::Vector3d> getLandmarks2D() {return m_landmarks2D;}
	std::vector<Eigen::Vector3d> getLandmarks3D() {return m_landmarks3D;}
	
private:
	void detect2DLandmarks()
	{
		string predictor_path="../template/shape_predictor_68_face_landmarks.dat";
		
		frontal_face_detector detector = get_frontal_face_detector();
		shape_predictor sp;
		deserialize(predictor_path) >> sp;
		
		dlib::cv_image<rgb_pixel> img(m_image);
		//load_image(img,image_path);
		
		//detect face
		std::vector<dlib::rectangle> dets = detector(img);
		dlib::rectangle det=dets[0];
	
		//detect landmarks
		full_object_detection shape=sp(img,det);
		

		int num = 68;
		m_landmarks2D.resize(num);
		for(size_t i=0; i<shape.num_parts();i++)
		{
			dlib::point p = shape.part(i);
			m_landmarks2D[i] = Eigen::Vector3d(p.x(),p.y(),1);

#ifdef VERBOSE
			cv::circle(m_image,cv::Point(p.x(),p.y()),20, cv::Scalar(0, 0, 255), -1, 4, 0);
#endif
		}
		
#ifdef VERBOSE
		cv::imwrite("2D landmark.jpg",m_image);
#endif
	}

	
	void detect3DLandmarks()
	{

		float threshold=50;  //m_image.rows/100;
		
		std::vector<trimesh::point> m_vertices = m_mesh->vertices;
		
		std::vector<std::vector<int> > projection(68);
		
		for(size_t i=0; i<m_vertices.size(); i++)
		{
			trimesh::point xyz=m_vertices[i];
			
			Eigen::Vector4d X3D(xyz[0],xyz[1],xyz[2],1);
			
			Eigen::Vector3d X2D;
			X2D = m_projMatrix*X3D;
			X2D /= X2D(2);
			
			for(int j=17; j<68; j++)
			{
				if(isNear(m_landmarks2D[j],X2D, threshold))
				{
					projection[j].push_back(i);
				}
			}
		}
		
#ifdef VERBOSE
		ofstream fout("3D landmark.obj");
#endif		
	
		for(int i=17; i<68; i++)  //17
		{
			int N = projection[i].size();
			std::vector<double> vecX,vecY,vecZ;
			
			
			//std::vector<std::pair<trimesh::point, int> > points;

			for(int j=0; j<N; j++)
			{
				int idx = projection[i][j];
				
				trimesh::point xyz=m_mesh->vertices[idx];
				
				vecX.push_back(xyz[0]);
				vecY.push_back(xyz[1]);
				vecZ.push_back(xyz[2]);
				
				//points.push_back(std::pair<trimesh::point, int>(xyz, idx));
			}
			
			//sort(points.begin(),points.end(), cmp);
			
			sort(vecX.begin(),vecX.end());
			sort(vecY.begin(),vecY.end());
			sort(vecZ.begin(),vecZ.end());

			m_landmarks3D.push_back(Eigen::Vector3d(vecX[N/2],vecY[N/2],vecZ[N/2]));
			
			//m_landmarks3D.push_back(Eigen::Vector3d(points[N/2].first[0],points[N/2].first[1],points[N/2].first[2]));
			//m_landmarks3D_Idx.push_back(points[N/2].second);
			
			
#ifdef VERBOSE		
			fout<<"v "<<vecX[N/2]<<" "<<vecY[N/2]<<" "<<vecZ[N/2]<<" 1 0 0"<<endl;
			//fout<<"v "<<points[N/2].first[0]<<" "<<points[N/2].first[0]<<" "<<points[N/2].first[0]<<" 1 0 0"<<endl;
#endif	
			
		}
		
#ifdef VERBOSE		
		fout.close();
#endif	
		
	}
	
};
