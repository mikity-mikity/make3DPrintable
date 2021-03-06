#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Random.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <ppl.h>
#include <array>
#include <vector>
#include <cassert>
#include <iostream>
#include <string>
#include <boost/math/constants/constants.hpp>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include<algorithm>	
#include<Eigen/Dense>
#include"GeometryProcessing.h"
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_3<K,CGAL::Fast_location> Delaunay;
typedef Delaunay::Point Point;
typedef CGAL::Vector_3<K> Vector;
typedef std::vector<Point> branch;
typedef struct
{
	Vector majorRadius;
	Vector minorRadius;
	Vector Normal;
}eclipse;
typedef struct
{
	Vector a;
}trapezoid;
using namespace std;
using namespace GeometryProcessing;
using namespace Concurrency;
typedef std::vector<eclipse> eclipses;
typedef boost::tuple<double,branch*> Rad_branch;
typedef boost::tuple<double,Mesh*> Thick_mesh;
typedef Delaunay::Cell_handle Cell_handle;
typedef Delaunay::Vertex_handle Vertex_handle;
typedef Delaunay::Locate_type    Locate_type;
double PI=0.0;

void computeGradient(Eigen::Vector3d &grad,vertex* v, int index, MeshStructure* MS, Eigen::VectorXd x, Eigen::Vector3d* normalPerFaces)
{
	grad(0)=0;
	grad(1)=0;
	grad(2)=0;
    for (auto e : v->onering)
    {
        auto N = normalPerFaces[e->owner->N];
        auto cx = x(index * 3 + 0);
        auto cy = x(index * 3 + 1);
        auto cz = x(index * 3 + 2);
        grad(0) += 2 * N.x() * N.x() * cx + 2 * N.x() * N.y() * cy + 2 * N.x() * N.z() * cz - 2 * N.x();
        grad(1) += 2 * N.y() * N.y() * cy + 2 * N.y() * N.z() * cz + 2 * N.y() * N.x() * cx - 2 * N.y();
        grad(2) += 2 * N.z() * N.z() * cz + 2 * N.z() * N.x() * cx + 2 * N.z() * N.y() * cy - 2 * N.z();
        double dot = cx * N.x() + cy * N.y() + cz * N.z();
		double norm = N.x()*N.x()+N.y()*N.y()+N.z()*N.z();
        grad(0) += 0.01*(2 * cx - 4 * dot * N.x() + norm * 2 * dot * N.x());
        grad(1) += 0.01*(2 * cy - 4 * dot * N.y() + norm * 2 * dot * N.y());
        grad(2) += 0.01*(2 * cz - 4 * dot * N.z() + norm * 2 * dot * N.z());
    }
}
void computeHessian(Eigen::Matrix3d &hess,vertex* v, int index,MeshStructure* MS, Eigen::VectorXd x, Vector3d* normalPerFaces)
{
	hess(0,0)=0;
	hess(0,1)=0;
	hess(0,2)=0;
	hess(1,0)=0;
	hess(1,1)=0;
	hess(1,2)=0;
	hess(2,0)=0;
	hess(2,1)=0;
	hess(2,2)=0;
    for (auto e : v->onering)
    {
        auto N = normalPerFaces[e->owner->N];
        auto cx = x(index * 3 + 0);
        auto cy = x(index * 3 + 1);
        auto cz = x(index * 3 + 2);
        hess(0,0)+= 2 * N.x() * N.x();
        hess(1,0)+= 2 * N.x() * N.y();
        hess(2,0)+= 2 * N.x() * N.z();
        hess(0,1)+= 2 * N.y() * N.x();
        hess(1,1)+= 2 * N.y() * N.y();
        hess(2,1)+= 2 * N.y() * N.z();
        hess(0,2)+= 2 * N.z() * N.x();
        hess(1,2)+= 2 * N.z() * N.y();
        hess(2,2)+= 2 * N.z() * N.z();
        double dot = cx * N.x() + cy * N.y() + cz * N.z();
		double norm = N.x()*N.x()+N.y()*N.y()+N.z()*N.z();
        hess(0,0) += 0.01*(2 - 4 * N.x() * N.x() + norm * 2 * N.x() * N.x());
        hess(1,0) += 0.01*(-4 * N.y() * N.x() + norm * 2 * N.y() * N.x());
        hess(2,0) += 0.01*(-4 * N.z() * N.x() + norm * 2 * N.z() * N.x());
        hess(1,1) += 0.01*(2 - 4 * N.y() * N.y() + norm * 2 * N.y() * N.y());
        hess(2,1) += 0.01*(-4 * N.z() * N.y() + norm * 2 * N.z() * N.y());
        hess(0,1) += 0.01*(-4 * N.x() * N.y() + norm * 2 * N.x() * N.y());
        hess(2,2) += 0.01*(2 - 4 * N.z() * N.z() + norm * 2 * N.z() * N.z());
        hess(0,2) += 0.01*(-4 * N.x() * N.z() + norm * 2 * N.x() * N.z());
        hess(1,2) += 0.01*(-4 * N.y() * N.z() + norm * 2 * N.y() * N.z());
    }
}

boost::tuple<double,GeometryProcessing::MeshStructure*,std::map<vertex*,boost::tuple<Eigen::Vector3d,Eigen::Vector3d>>> computeNormal(boost::tuple<double,Mesh*,GeometryProcessing::MeshStructure* >input)
{
	Mesh *myMesh;
	GeometryProcessing::MeshStructure *MS;
	double t;
	boost::tie(t,myMesh,MS)=input;
    __int64 numvar = MS->vertices.size()* 3;
	
    Eigen::VectorXd x = Eigen::MatrixXd::Zero(numvar, 1);
    vector<Point3d> nodes = myMesh->Vertices;
    Eigen::Vector3d *normalPerFaces = new Eigen::Vector3d[MS->nFaces()];
    //initial guess
    for (auto v : MS->vertices)
    {
		std::vector<vertex*>::iterator iter = std::find(MS->vertices.begin(), MS->vertices.end(), v);
		__int64 index = std::distance(MS->vertices.begin(), iter);
		
        Eigen::Vector3d N=Eigen::Vector3d(0,0,0);
        for (auto e : v->onering)
        {
            //compute normal
            auto P = nodes[e->P->N];
            auto Q = nodes[e->next->P->N];
            auto R = nodes[e->next->next->P->N];
            auto b = Eigen::Vector3d(P.X - Q.X,P.Y-Q.Y,P.Z-Q.Z);
            auto c = Eigen::Vector3d(R.X - Q.X,R.Y-Q.Y,R.Z-Q.Z);
            auto n = b.cross(c);
			n.normalize();
            normalPerFaces[e->owner->N] = n;
            N += n;
        }
        N.normalize();
        x(index * 3 + 0) = N(0);
        x(index * 3 + 1) = N(1);
        x(index * 3 + 2) = N(2);
    }
	double tol=0.00000001;
	int index=0;
	Eigen::Vector3d grad=Eigen::Vector3d(0,0,0);
	Eigen::Matrix3d hess=Eigen::Matrix3d().Zero();
	for(auto itr=MS->vertices.begin();itr!=MS->vertices.end();itr++,index++)
	{
		for (int i = 0; i < 500; i++)
		{
			computeGradient(grad,*itr,index, MS,x, normalPerFaces);
			computeHessian(hess,*itr,index,MS, x, normalPerFaces);
			Eigen::PartialPivLU<Matrix3d> sol(hess);
			Eigen::Vector3d dx = sol.solve(-grad);
			x(index * 3 + 0) += dx(0);
            x(index * 3 + 1) += dx(1);
            x(index * 3 + 2) += dx(2);

			if (grad.norm()< tol)break;
		}
	}
	std::map<vertex*,boost::tuple<Eigen::Vector3d,Eigen::Vector3d>> output;
	index=0;
	for(auto itr=MS->vertices.begin();itr!=MS->vertices.end();itr++,index++)
	{
		auto v=*itr;	
		Eigen::Vector3d P(nodes[v->N].X,nodes[v->N].Y,nodes[v->N].Z);
        Eigen::Vector3d N(x(index * 3 + 0),x(index * 3 + 1),x(index * 3 + 2));
        output.insert(std::pair<vertex*,boost::tuple<Eigen::Vector3d,Eigen::Vector3d>>(v,boost::make_tuple(P,N)));
	}
	delete[] normalPerFaces;
	return boost::make_tuple(t,MS,output);
}
boost::tuple<double,double> read(string filename,std::vector<Rad_branch> &data,std::vector<Thick_mesh> &mData)
{
	branch* _branch;
	Mesh* _mesh;
	string line;
	double R;
	double T;
	int N;
	int nV,nF;
	// read file
	double minR=10000;
	double minT=10000;

	ifstream ifs(filename);  //input file
	if(ifs.fail()){
		std::cout<<"File does not exist."<<endl;
		exit(0);
	}

	while(getline(ifs,line))
	{
		std::vector<string> words;
		boost::algorithm::split( words, line, boost::algorithm::is_any_of(","));
		char _prefix[20];
		sscanf(words[0].data(),"%s",_prefix);
		string prefix=string(_prefix);
		if(prefix=="P")
		{
			_branch=new branch();
			sscanf(words[1].data(),"%d",&N);
		}
		if(prefix=="R")
		{
			sscanf(words[1].data(),"%lf",&R);
			if(R<minR)minR=R;
			data.push_back(boost::make_tuple(R,_branch));
		}
		if(prefix=="C"){
			double x,y,z;
			sscanf(words[1].data(),"%lf",&x);
			sscanf(words[2].data(),"%lf",&y);
			sscanf(words[3].data(),"%lf",&z);
			_branch->push_back(Point(x,y,z));
		}
		if(prefix=="S")
		{
			_mesh=new Mesh();
			sscanf(words[1].data(),"%d",&nV);
			sscanf(words[2].data(),"%d",&nF);
		}
		if(prefix=="T")
		{
			sscanf(words[1].data(),"%lf",&T);
			if(T<minT)minT=T;
			mData.push_back(boost::make_tuple(T,_mesh));
		}
		if(prefix=="V"){
			double x,y,z;
			sscanf(words[1].data(),"%lf",&x);
			sscanf(words[2].data(),"%lf",&y);
			sscanf(words[3].data(),"%lf",&z);
			_mesh->Vertices.push_back(Point3d(x,y,z));
		}
		if(prefix=="F"){
			int A,B,C;
			sscanf(words[1].data(),"%d",&A);
			sscanf(words[2].data(),"%d",&B);
			sscanf(words[3].data(),"%d",&C);
			_mesh->Faces.push_back(MeshFace(A,B,C));
		}
		
	}
	ifs.close();
	return boost::make_tuple(minR,minT);
}
void generate_eclipseTree(std::vector<Rad_branch> &data,std::vector<eclipses*> &eclipseTree)
{
	int C=0;
	int num=0;

	//Construct eclipses at each node	
	for(vector<Rad_branch>::iterator itr=data.begin();itr!=data.end();itr++,C++)
	{
		double Radius;
		branch* _branch;
		boost::tie(Radius,_branch)=*itr;
		eclipses *_eclipses=new eclipses();
		eclipseTree.push_back(_eclipses);
		for(branch::iterator anitr=_branch->begin();anitr!=_branch->end();anitr++)
		{
			eclipse newEclipse;
			if(anitr==_branch->begin())
			{
				Point P=*(anitr);
				Point Q=*(anitr+1);
				Vector V=(Q-P);
				V=V/std::sqrt(V.squared_length());
				Vector Z(0,0,1);
				double t=V*Z;
				if(std::fabs(t)>0.9)
				{
					Z=Vector(0,1,0);
				}
				Vector W=CGAL::cross_product(V,Z);
				Vector T=CGAL::cross_product(W,V);
				W=W/std::sqrt(W.squared_length());
				T=T/std::sqrt(T.squared_length());
				newEclipse.majorRadius=W;
				newEclipse.minorRadius=T;
				newEclipse.Normal=CGAL::cross_product(W,T);
			}else if(anitr==_branch->end()-1)
			{
				Point Q=*(anitr-1);
				Point R=*(anitr);
				Vector V=(R-Q);
				V=V/std::sqrt(V.squared_length());
				Vector Z(0,0,1);
				double t=V*Z;
				if(std::fabs(t)>0.9)
				{
					Z=Vector(0,1,0);
				}
				Vector W=CGAL::cross_product(V,Z);
				Vector T=CGAL::cross_product(W,V);
				W=W/std::sqrt(W.squared_length());
				T=T/std::sqrt(T.squared_length());
				newEclipse.majorRadius=W;
				newEclipse.minorRadius=T;
				newEclipse.Normal=CGAL::cross_product(W,T);
			}else
			{
				Point P=*(anitr-1);
				Point Q=*(anitr);
				Point R=*(anitr+1);
				Vector V1=(Q-P);
				Vector V2=(R-Q);
				V1=V1/std::sqrt(V1.squared_length());
				V2=V2/std::sqrt(V2.squared_length());
				if(V1*V2==1.0)//if parallel
				{
					Vector Z(0,0,1);
					double t=V1*Z;
					if(std::fabs(t)>0.9)
					{
						Z=Vector(0,1,0);
					}
					Vector W=CGAL::cross_product(V1,Z);
					Vector T=CGAL::cross_product(W,V1);
					W=W/std::sqrt(W.squared_length());
					T=T/std::sqrt(T.squared_length());
					newEclipse.majorRadius=W;
					newEclipse.minorRadius=T;
					newEclipse.Normal=CGAL::cross_product(W,T);
				}else
				{
					//angle between two lines
					double inner=V1*V2;
					double alpha=std::acos(inner);
					double theta=alpha/2.;

					Vector N=CGAL::cross_product(V1,V2);
					Vector E=CGAL::cross_product(V2,N);
					E=E/std::sqrt(E.squared_length());
					Vector R=(std::cos(theta)*E-std::sin(theta)*V2);
					newEclipse.majorRadius=R;
					newEclipse.minorRadius=N;
					newEclipse.Normal=CGAL::cross_product(R,N);
				}

			}
			_eclipses->push_back(newEclipse);
		}
	}
}
inline double max(std::vector<double> v)
{
    double maxVal = -1000000;       // �����ŏ��l
    for(int i = 0; i < (int)v.size(); ++i) {
        if( v[i] > maxVal )
            maxVal = v[i];
    }
    return maxVal;
}
void generate_exterior2(boost::tuple<double,GeometryProcessing::MeshStructure*,std::map<vertex*,boost::tuple<Eigen::Vector3d,Eigen::Vector3d>>> tMS,double baseRes,std::vector<Point> &exterior)
{
	double t;
	GeometryProcessing::MeshStructure *MS;
	std::map<vertex*,boost::tuple<Eigen::Vector3d,Eigen::Vector3d>> info;
	boost::tie(t,MS,info)=tMS;
	double tt[5]={-0.5,-0.25,0.0,0.25,0.5};
	for(auto v:MS->vertices)
	{
		if(!v->isBoundary()){
			Eigen::Vector3d P,N; //Position,Normal
			boost::tie(P,N)=info.find(v)->second;
			for(int i=0;i<5;i++)
			{
				Eigen::Vector3d Pc=P+tt[i]*t*N;
				exterior.push_back(Point(Pc.x(),Pc.y(),Pc.z()));
			}
		}
	}
	for(auto e:MS->edges())
	{
		if((!e->isBoundary())&&(!e->ifPairIsBoundary()))
		{
			Eigen::Vector3d P1,P2,N1,N2; //Position,Normal
			boost::tie(P1,N1)=info.find(e->P)->second;
			boost::tie(P2,N2)=info.find(e->next->P)->second;
			for(int i=0;i<5;i++)
			{

				Eigen::Vector3d P1in=P1+tt[i]*t*N1;
				Eigen::Vector3d P2in=P2+tt[i]*t*N2;
				auto T=P2in-P1in;
				double L=T.norm();
				int LDIV=2;
				if(L>baseRes*2.0)LDIV=(int)(L/baseRes);
				for(int ss=1;ss<LDIV;ss++)
				{
					double sss=((double)ss)/((double)LDIV);
					Eigen::Vector3d Pc=P1in+T*sss;
					exterior.push_back(Point(Pc.x(),Pc.y(),Pc.z()));
				}
			}
		}
	}
	for(auto f:MS->faces)
	{
		Eigen::Vector3d P1,P2,P3,N1,N2,N3,P1out,P1in,P2out,P2in,P3out,P3in,_P1,_P2,_P3; //Position,Normal
		boost::tie(P1,N1)=info.find(f->firsthalfedge->P)->second;
		boost::tie(P2,N2)=info.find(f->firsthalfedge->next->P)->second;
		boost::tie(P3,N3)=info.find(f->firsthalfedge->next->next->P)->second;
		P1in=P1+N1*t*0.5;
		P2in=P2+N2*t*0.5;
		P3in=P3+N3*t*0.5;
		P1out=P1-N1*t*0.5;
		P2out=P2-N2*t*0.5;
		P3out=P3-N3*t*0.5;
		auto T12in=P2in-P1in;
		auto T23in=P3in-P2in;
		auto T31in=P1in-P3in;
		auto T12out=P2out-P1out;
		auto T23out=P3out-P2out;
		auto T31out=P1out-P3out;
		//choose longest
		double dd[]={T12in.norm(),T23in.norm(),T31in.norm(),T12out.norm(),T23out.norm(),T31out.norm()};
		std::vector<double> ll(dd,dd+6);
		double L=max(ll);
		int LDIV=(int)(L/baseRes);
		if(LDIV<2)LDIV=2;
		for(int i=0;i<5;i++)
		{
			_P1=P1+N1*t*tt[i];
			_P2=P2+N2*t*tt[i];
			_P3=P3+N3*t*tt[i];
			//choose longest
			//auto T12=_P2-_P1;
			auto T23=_P3-_P2;
			auto T13=_P3-_P1;
			for(int v=1;v<LDIV;v++)
			{
				double sv=((double)v)/((double)LDIV);
				auto __P1=_P1+T13*sv;
				auto __P2=_P2+T23*sv;
				for(int u=1;u<LDIV-v;u++)
				{
					double su=((double)u)/((double)(LDIV-v));
					auto T12=__P2-__P1;
					auto p=__P1+T12*su;
					exterior.push_back(Point(p.x(),p.y(),p.z()));
				}
			}
		}
	}
	
	auto e=MS->boundaryStart;
	do{
		if(!e->isBoundary())std::cout<<"error!"<<endl;
		Eigen::Vector3d P,Q,Pin,Qin,Pout,Qout,Pt,Qt,PQ,p,N1,N2; //Position,Normal
		boost::tie(P,N1)=info.find(e->P)->second;
		boost::tie(Q,N2)=info.find(e->next->P)->second;
		Pin=P+N1*t*0.5;
		Qin=Q+N2*t*0.5;
		Pout=P-N1*t*0.5;
		Qout=Q-N2*t*0.5;
		//choose longest
		Eigen::Vector3d PQin=Qin-Pin;
		Eigen::Vector3d PQout=Qout-Pout;
		double uL=std::max(PQin.norm(),PQout.norm());
		double vL=std::max(N1.norm(),N2.norm());
		int UDIV=(int)(uL/baseRes);
		int VDIV=(int)(vL/baseRes);
		if(UDIV<2)UDIV=2;
		if(VDIV<4)VDIV=4;
		for(int vs=0;vs<=VDIV;vs++)
		{
			double vss=((double)vs)/((double)VDIV);
			Pt=P+N1*t*(vss-0.5);
			Qt=Q+N2*t*(vss-0.5);
			PQ=Qt-Pt;
			for(int us=0;us<UDIV;us++)
			{
				double uss=((double)us)/((double)UDIV);
				p=Pt+PQ*uss;
				exterior.push_back(Point(p.x(),p.y(),p.z()));
			}
		}
		e=e->next;
	}while(e!=MS->boundaryStart);
	
}
void generate_exterior(std::vector<Rad_branch> &data,std::vector<eclipses*> &eclipseTree,double baseRes,std::vector<Point> &exterior)
{
	vector<Rad_branch>::iterator itrA=data.begin();
	vector<eclipses*>::iterator itrB=eclipseTree.begin();

	__int64 NN=(data.size()/20);
	if(NN<1)NN=1;	
	int N=0;
	while(itrA!=data.end())
	{
		double Radius;
		branch* _branch;
		eclipses* _eclipses=*itrB;
		boost::tie(Radius,_branch)=*itrA;
		branch::iterator itrC=_branch->begin();
		eclipses::iterator itrD=_eclipses->begin();
		Vector beforeX(0,0,0);
		Vector beforeY(0,0,0);
		//Division number along diameter
		int RDIV=12;
		if(Radius*2.*PI/12.<baseRes)
		{
			RDIV=12;
		}else{
			RDIV=(int)(Radius*2.*PI/baseRes);
		}
		while(itrC!=_branch->end()-1)
		{
			//Division number along line
			Point P=*itrC;
			Point Q=*(itrC+1);
			Vector V=Q-P;
			double Length=std::sqrt(V.squared_length());
			V=V/Length;
			int DIV=1;//default division number
			if(Length<baseRes)
			{
				DIV=1;
			}else{
				DIV=(int)(Length/baseRes)+1;
			}
			//if in the first run, before is igonored
			//from the second run, before is used

			if(itrC==_branch->begin())
			{
				Vector Z(0,0,1);
				double t=V*Z;
				if(std::fabs(t)>0.9)
				{
					Z=Vector(0,1,0);
				}
				beforeX=CGAL::cross_product(V,Z);
				beforeY=CGAL::cross_product(V,beforeX);
				beforeX=beforeX/std::sqrt(beforeX.squared_length());
				beforeY=beforeY/std::sqrt(beforeY.squared_length());
			}else
			{				
				//project beforeX and beforeY to the plane at the node
				double cx=-(beforeX*itrD->Normal)/(V*itrD->Normal);
				double cy=-(beforeY*itrD->Normal)/(V*itrD->Normal);
				Vector tmpX=beforeX+cx*V;
				Vector tmpY=beforeY+cy*V;

				//Compute plane	
				Vector Z(0,0,1);
				double t=V*Z;
				if(std::fabs(t)>0.9)
				{
					Z=Vector(0,1,0);
				}
				Vector X=CGAL::cross_product(V,Z);
				Vector Y=CGAL::cross_product(V,beforeX);
				Vector N=CGAL::cross_product(X,Y);
				//project tmpX and tmpY to the plane conformed of X and Y
				beforeX=tmpX-(tmpX*N)*N;
				beforeY=CGAL::cross_product(V,beforeX);
				beforeX=beforeX/std::sqrt(beforeX.squared_length());
				beforeY=beforeY/std::sqrt(beforeY.squared_length());
			}
			int DIV2=DIV;
			if(itrC==_branch->end()-2)DIV2=DIV+1;
			for(int ss=0;ss<DIV2;ss++)
			{
				double s=((double)ss)/((double)DIV);
				exterior.push_back(Point((*itrC).x()*(1-s)+(*(itrC+1)).x()*s,(*itrC).y()*(1-s)+(*(itrC+1)).y()*s,(*itrC).z()*(1-s)+(*(itrC+1)).z()*s));
				for(int i=0;i<RDIV;i++)
				{
					double theta=(double)(i)/RDIV*2.*PI;

					Vector BI=0.8*Radius*(beforeX*std::cos(theta)+beforeY*std::sin(theta));
					Vector BE=1.05*Radius*(beforeX*std::cos(theta)+beforeY*std::sin(theta));

					//Project N
					double cI1=-(BI*itrD->Normal)/(V*itrD->Normal);
					double cI2=-(BI*(itrD+1)->Normal)/(V*(itrD+1)->Normal);
					double cE1=-(BE*itrD->Normal)/(V*itrD->Normal);
					double cE2=-(BE*(itrD+1)->Normal)/(V*(itrD+1)->Normal);
					Vector tmpI1=BI+cI1*V;
					Vector tmpI2=BI+cI2*V;
					Vector tmpE1=BE+cE1*V;
					Vector tmpE2=BE+cE2*V;
					Point DI1=(*itrC)+tmpI1;
					Point DI2=(*(itrC+1))+tmpI2;
					Point DE1=(*itrC)+tmpE1;
					Point DE2=(*(itrC+1))+tmpE2;
					if(itrC==_branch->end()-2)
					{
						Point DI(DI2.x()*s+DI1.x()*(1-s),DI2.y()*s+DI1.y()*(1-s),DI2.z()*s+DI1.z()*(1-s));
						Point DE(DE2.x()*s+DE1.x()*(1-s),DE2.y()*s+DE1.y()*(1-s),DE2.z()*s+DE1.z()*(1-s));
						exterior.push_back(DI);
						exterior.push_back(DE);
					}else
					{
						Point DI(DI2.x()*s+DI1.x()*(1-s),DI2.y()*s+DI1.y()*(1-s),DI2.z()*s+DI1.z()*(1-s));
						Point DE(DE2.x()*s+DE1.x()*(1-s),DE2.y()*s+DE1.y()*(1-s),DE2.z()*s+DE1.z()*(1-s));
						exterior.push_back(DI);
						exterior.push_back(DE);
					}
				}
			}
			itrC++;
			itrD++;
		}

		itrA++;
		itrB++;
		if(((int)N/NN)*NN==N)
		{
			std::cout<<"*";
		}
		N++;
	}
	std::cout<<endl;
}
void triangulate(Delaunay &T,std::vector<Point> &exterior)
{
	for(int i=0;i<20;i++)
	{
		std::vector<Point>::iterator be=exterior.begin()+exterior.size()*i/20;
		std::vector<Point>::iterator en=exterior.begin()+exterior.size()*(i+1)/20;
		T.insert(be,en);
		std::cout<<"*";
	}
	exterior.clear();
	std::cout<<endl;

}
void asign_index_to_cells(std::map<Delaunay::Cell_handle,__int64> &index,Delaunay &T)
{
	__int64 N=0;
	for(Delaunay::Finite_cells_iterator itr=T.finite_cells_begin();itr!=T.finite_cells_end();itr++)
	{
		index.insert(pair<const Delaunay::Finite_cells_iterator,__int64>(itr,N));
		N++;
	}
}
void init_cells(Delaunay &T,__int64* cells)
{
	__int64 S=T.number_of_finite_cells();
	__int64* ptr1=cells;
	for(int i=0;i<S;i++)
	{
		*ptr1=0;
		ptr1++;
	}
}
typedef struct
{
	int DIV;
	Vector beforeX;
	Vector beforeY;
	Vector V;
	Cell_handle before;
	vector<Point> particles;
	double s;
	branch::iterator itrC;
	eclipses::iterator itrD;
	branch* branch;
}info;
typedef struct
{
	Cell_handle before;
	Point p;
}info2;


__int64 computeInterior2(vector<boost::tuple<double,GeometryProcessing::MeshStructure*,std::map<vertex*,boost::tuple<Eigen::Vector3d,Eigen::Vector3d>>>> mesh_infos,double baseRes,std::function<void(info2&,__int64&)>func,__int64 nsize)
{	
	__int64 numInterior=0;
	__int64 next=0;
	double tt[5]={-0.5,-0.25,0.0,0.25,0.5};
	double uv_ser[] = { 0.02, 0.05, 0.1, 0.22, 0.39, 0.61, 0.78, 0.9, 0.95, 0.98 };
	double t_ser[40];
	for(int i=0;i<4;i++)
	{
		for(int j=0;j<10;j++)
		{
			int index=i*10+j;
			double val=tt[i]+(tt[i+1]-tt[i])*(((double)j+0.5)/10);
			t_ser[index]=val;
		}
	}
	for (auto tMS : mesh_infos)
	{
		double t;
		GeometryProcessing::MeshStructure *MS;
		std::map<vertex*,boost::tuple<Eigen::Vector3d,Eigen::Vector3d>> info;
		boost::tie(t,MS,info)=tMS;
		int TDIV=40;
		vector<boost::tuple<vector<face*>::iterator, vector<face*>::iterator, info2*, int>> tasks;
		__int64 size= MS->faces.size();
		int nTasks = 100;
		if (size < 1000)nTasks = 1;
		info2* myInfo=new info2[nTasks];
		
		for (int i = 0; i < nTasks; i++)
		{
			__int64 be = (size*i) / nTasks;
			__int64 en = (size*(i + 1)) / nTasks;
			myInfo[i].before = NULL;
			vector<face*>::iterator first = MS->faces.begin(), second = MS->faces.begin();
			for (int j = 0; j < be; j++)
			{
				first++;
			}
			for (int j = 0; j < en; j++)
			{
				second++;
			}
			tasks.push_back(boost::make_tuple<vector<face*>::iterator, vector<face*>::iterator, info2*,int>(first, second, &myInfo[i], i));
		}
		critical_section cs;
		parallel_for_each(tasks.begin(), tasks.end(), [&baseRes, &TDIV, &t, &func, &info, &numInterior, &next, &nsize,&uv_ser,&t_ser,&cs](boost::tuple<vector<face*>::iterator, vector<face*>::iterator, info2*, int> tup)
		{
			vector<face*>::iterator begin;
			vector<face*>::iterator end;
			info2* myInfo;
			int taskNum;
			boost::tie(begin, end, myInfo, taskNum) = tup;
			for (vector<face*>::iterator itr = begin; itr != end; itr++)
			{
				face *f = *itr;
				Eigen::Vector3d P1, P2, P3, N1, N2, N3, P1out, P1in, P2out, P2in, P3out, P3in, _P1, _P2, _P3; //Position,Normal
				boost::tie(P1, N1) = info.find(f->firsthalfedge->P)->second;
				boost::tie(P2, N2) = info.find(f->firsthalfedge->next->P)->second;
				boost::tie(P3, N3) = info.find(f->firsthalfedge->next->next->P)->second;
				P1in = P1 + N1*t*0.5;
				P2in = P2 + N2*t*0.5;
				P3in = P3 + N3*t*0.5;
				P1out = P1 - N1*t*0.5;
				P2out = P2 - N2*t*0.5;
				P3out = P3 - N3*t*0.5;
				auto T12in = P2in - P1in;
				auto T23in = P3in - P2in;
				auto T31in = P1in - P3in;
				auto T12out = P2out - P1out;
				auto T23out = P3out - P2out;
				auto T31out = P1out - P3out;
				//choose longest
				double dd[] = { T12in.norm(), T23in.norm(), T31in.norm(), T12out.norm(), T23out.norm(), T31out.norm() };
				std::vector<double> ll(dd, dd + 6);
				double L = max(ll);
				int LDIV = (int)(L / baseRes);
				if (LDIV < 2)LDIV = 2;
				//LDIV *= 10;
				__int64 nI = 0;
				for (int ss = 0; ss <TDIV; ss++)
				{
					double sss = t_ser[ss];
					_P1 = P1 + N1*t*sss;
					_P2 = P2 + N2*t*sss;
					_P3 = P3 + N3*t*sss;
					//choose longest
					//auto T12=_P2-_P1;
					auto T23 = _P3 - _P2;
					auto T13 = _P3 - _P1;
					for (int v = 0; v < LDIV; v++)
					{
						for (int vv = 0; vv < 10; vv++)
						{
							double sv = ((double)v /*- 0.5*/ + uv_ser[vv]) / ((double)LDIV);
							auto __P1 = _P1 + T13*sv;
							auto __P2 = _P2 + T23*sv;
							for (int u = 0; u < LDIV - v; u++)
							{
								for (int uu = 0; uu < 10; uu++)
								{
									double su = ((double)u /*- 0.5*/ + uv_ser[uu]) / ((double)(LDIV - v));
									auto T12 = __P2 - __P1;
									auto p = __P1 + T12*su;
									myInfo->p = Point(p.x(), p.y(), p.z());
									func(*myInfo, nI);
								}
							}
						}
					}
				}
				cs.lock();
				numInterior += nI;
				while (numInterior >= next)
				{
					next += 10000000;
					if (nsize == 0)
						std::cout << taskNum << ":" << next - 10000000 << endl;
					else
						std::cout << taskNum << ":" << next - 10000000 << "/" << nsize << endl;

				}
				cs.unlock();
			}
		});
		delete [] myInfo;
	}
	//std::cout << "numInterior2:" << numInterior << endl;
	//std::cout << "next:" << next << endl;
	//std::cin.get();
	return numInterior;
}
__int64 computeInterior(std::vector<Rad_branch> &data, std::vector<eclipses*> &eclipseTree, double baseRes, std::function<void(info&,__int64&)>func,__int64 nsize)
{

	vector<Rad_branch>::iterator __itrA = data.begin();
	vector<eclipses*>::iterator __itrB = eclipseTree.begin();
	__int64 numInterior=0;
	__int64 next=0;
	info* __myInfo=new info[data.size()];
	int N = 0;
	vector<boost::tuple<std::vector<Rad_branch>::iterator,std::vector<eclipses*>::iterator, info*, int>> tasks;
	for (int i = 0; i < data.size(); i++)
	{
		__myInfo[i].before = NULL;
		if (__itrA == data.end())break;
		if (__itrB == eclipseTree.end())break;
		tasks.push_back(boost::make_tuple(__itrA, __itrB, &__myInfo[i], i));
		__itrA++;
		__itrB++;
	}
	critical_section cs;

    parallel_for_each(tasks.begin(), tasks.end(), [&baseRes, &numInterior, &next, &func,&nsize,&cs](boost::tuple<std::vector<Rad_branch>::iterator, std::vector<eclipses*>::iterator, info*, int> task)
		//while(itrA!=data.end())
	{
		int YDIV = 8;
		double Radius;
		std::vector<Rad_branch>::iterator itrA;
		std::vector<eclipses*>::iterator itrB;
		info* myInfo;
		int taskNum;
		boost::tie(itrA, itrB, myInfo, taskNum) = task;
		eclipses* _eclipses = *itrB;
		boost::tie(Radius, myInfo->branch) = *itrA;  //decompose
		myInfo->itrC = myInfo->branch->begin();
		myInfo->itrD = _eclipses->begin();
		myInfo->beforeX = Vector(0, 0, 0);
		myInfo->beforeY = Vector(0, 0, 0);
		//Division number along diameter
		int RDIV = 12;
		if (Radius*2.*PI / 12. < baseRes)
		{
			RDIV = 12;
		}
		else{
			RDIV = (int)(Radius*2.*PI / baseRes);
		}

		//Generate base particles
		vector<Vector> vectors;
		double alpha = PI / ((double)RDIV);
		double R2 = Radius*std::cos(alpha);
		for (double i = 0.5; i < RDIV; i++)
		{
			double theta = 2.*PI*((double)i) / ((double)RDIV);
			Vector vector(R2*std::cos(theta), R2*std::sin(theta), 0);
			vectors.push_back(vector);
		}
		double edge = 2 * std::sin(alpha)*Radius;
		double dx = edge / YDIV;
		myInfo->particles.clear();

		for (double i = -Radius; i < Radius; i += dx)
		{
			for (double j = -Radius; j < Radius; j += dx)
			{
				Vector V(i, j, 0);
				bool flag = true;
				for (auto T = vectors.begin(); T != vectors.end(); T++)
				{
					if ((*T)*V < R2*R2)
					{
						continue;
					}
					else
					{
						flag = false;
						break;
					}

				}
				if (flag)
					myInfo->particles.push_back(Point(i, j, 0));
			}
		}

		std::cout << taskNum << ":particles.size()" << myInfo->particles.size() << endl;
		std::cout << taskNum << ":_branch.size()" << myInfo->branch->size() << endl;

		while (myInfo->itrC != myInfo->branch->end() - 1)
		{
			//Division number along line
			Point P = *myInfo->itrC;
			Point Q = *(myInfo->itrC+1);
			myInfo->V = Q - P;
			double Length = std::sqrt(myInfo->V.squared_length());
			myInfo->V = myInfo->V / Length;
			myInfo->DIV = 1;//default division number
			if (Length < baseRes)
			{
				myInfo->DIV = 1;
			}
			else{
				myInfo->DIV = (int)(Length / baseRes) + 1;
			}
			myInfo->DIV *= 5;
			__int64 __nI = 0;
			func(*myInfo,__nI);
			cs.lock();
			numInterior += __nI;
			while(numInterior >= next)
			{
				next += 10000000;
				if (nsize==0)
					std::cout << taskNum << ":" << next - 10000000 << endl;
				else
					std::cout << taskNum << ":" << next - 10000000 << "/" << nsize << endl;
			}
			cs.unlock();
			myInfo->itrC++;
			myInfo->itrD++;
		}
	});
	delete [] __myInfo;
	//std::cout << "numInterior1:" << numInterior << endl;
	//std::cout << "next:" << next << endl;
	//std::cin.get();
	return numInterior;
}
int main(int argc, char *argv[])
{
	if(argc<2)return 0;
	string filename=argv[1];   //input filename
    std::vector<string> left_right;
	
	boost::algorithm::split(left_right , filename, boost::algorithm::is_any_of("."));
	string NAME=left_right[0];

	PI=boost::math::constants::pi<double>();
	std::cout<<"start reading file"<<"["<<filename<<"]"<<endl;


	//File write	
	string filename1=NAME+".out";    //vertices
	string filename2=NAME+"_S.out";  //indices of coner viertices of triangles
	string filename3=NAME+"_F.out";  //tets
	//string filename4=NAME+"_D.out";  //interior points
	ofstream ofs(filename1);
	ofstream ofs2(filename2);
	ofstream ofs3(filename3);
	//ofstream ofs4(filename4);	

	std::vector<Rad_branch> data;
	std::vector<Thick_mesh> mData;
	double minR,minT;
	boost::tie(minR,minT)=read(filename,data,mData);
	double baseRes=std::min(minT/4.,minR*2*PI/12.);
	cout<<"baseRes="<<baseRes<<endl;
	std::vector<boost::tuple<double,Mesh*,GeometryProcessing::MeshStructure*>> meshStructures;

	for(auto tM:mData)
	{
		double t;
		Mesh *m;
		boost::tie(t,m)=tM;
		GeometryProcessing::MeshStructure *MS=GeometryProcessing::MeshStructure::CreateFrom(m);
		meshStructures.push_back(boost::make_tuple(t,m,MS));
	}
	vector<boost::tuple<double,GeometryProcessing::MeshStructure*,std::map<vertex*,boost::tuple<Eigen::Vector3d,Eigen::Vector3d>>>> mesh_infos;
	for(auto MS:meshStructures)
	{
		mesh_infos.push_back(computeNormal(MS));
	}
	std::vector<Point> exterior;
	std::vector<eclipses*> eclipseTree;
	generate_eclipseTree(data,eclipseTree);

	std::cout<<"construct exterior"<<endl;
	for(auto tMS:mesh_infos)
	{
		generate_exterior2(tMS,baseRes,exterior);
	}
	generate_exterior(data,eclipseTree,baseRes,exterior);
	std::cout<<"exterior.size()"<<exterior.size()<<endl;

	
	// building their Delaunay triangulation.
	std::cout<<"start triangulation"<<endl;
	Delaunay T;
	triangulate(T,exterior);
	std::cout<<"end triangulation"<<endl;
	exterior.clear();
	std::cout<<"T.number_of_vertices:"<<T.number_of_vertices()<<endl;
	std::cout<<"number_of_finite_cells:"<<T.number_of_finite_cells()<<endl;
	std::map<Delaunay::Cell_handle,__int64> index;
	std::cout<<"create index"<<endl;
	asign_index_to_cells(index,T);
	std::cout<<"start cell search"<<endl;

	
	__int64* cells=new __int64[T.number_of_finite_cells()];
	
	std::cout<<"initialize cells"<<endl;
	init_cells(T,cells);
	std::cout << "1" << endl;

	
	
	//compute total number of interior points

	__int64 size=computeInterior(data,eclipseTree,baseRes,[](info& myInfo,__int64& numInterior){
		for(double ss=0.5;ss<myInfo.DIV;ss++)
			{
				for(auto particle=myInfo.particles.begin();particle!=myInfo.particles.end();particle++)
				{
					++numInterior;
				}
			}
	}, 0);
	std::cout<<"interior.size():"<< size<<endl;
	__int64 size2 = computeInterior2(mesh_infos, baseRes, [](info2& myInfo, __int64& numInterior){
		++	numInterior;
	},0);
	std::cout << "interior.size():" << size << endl;
	std::cout << "interior2.size():" << size2 << endl;
	computeInterior2(mesh_infos, baseRes, [&cells, &index, &T](info2& myInfo, __int64& numInterior){
		int li, lj;
		Locate_type lt;
		Cell_handle c;
		if (myInfo.before==NULL)
			c = T.locate(myInfo.p, lt,li,lj);
		else
			c=T.locate(myInfo.p,lt,li,lj,myInfo.before);
		if(lt==Locate_type::CELL)
		{
				std::map<Delaunay::Cell_handle,__int64>::iterator it_c=index.find(c);
				__int64 d=it_c->second;
				cells[d]++;
				myInfo.before=c;
		}
					
		++numInterior;
	}, size2);
	for(auto tMS:meshStructures)
	{
		double t;
		Mesh *mesh;
		MeshStructure *MS;
		boost::tie(t,mesh,MS)=tMS;
		delete MS;
		delete mesh;
	}
	meshStructures.clear();
	computeInterior(data, eclipseTree, baseRes, [&cells, &index, &T](info& myInfo, __int64& numInterior){

			Cell_handle c;
					
			if(myInfo.itrC==myInfo.branch->begin())
			{
				Vector Z(0,0,1);
				double t=myInfo.V*Z;
				if(std::fabs(t)>0.9)
				{
					Z=Vector(0,1,0);
				}
				myInfo.beforeX=CGAL::cross_product(myInfo.V,Z);
				myInfo.beforeY=CGAL::cross_product(myInfo.V,myInfo.beforeX);
				myInfo.beforeX=myInfo.beforeX/std::sqrt(myInfo.beforeX.squared_length());
				myInfo.beforeY=myInfo.beforeY/std::sqrt(myInfo.beforeY.squared_length());
			}else
			{				
				//project beforeX and beforeY to the plane at the node
				double cx=-(myInfo.beforeX*myInfo.itrD->Normal)/(myInfo.V*myInfo.itrD->Normal);
				double cy=-(myInfo.beforeY*myInfo.itrD->Normal)/(myInfo.V*myInfo.itrD->Normal);
				Vector tmpX=myInfo.beforeX+cx*myInfo.V;
				Vector tmpY=myInfo.beforeY+cy*myInfo.V;

				//Compute plane	
				Vector Z(0,0,1);
				double t=myInfo.V*Z;
				if(std::fabs(t)>0.9)
				{
					Z=Vector(0,1,0);
				}
				Vector X=CGAL::cross_product(myInfo.V,Z);
				Vector Y=CGAL::cross_product(myInfo.V,myInfo.beforeX);
				Vector N=CGAL::cross_product(X,Y);
				//project tmpX and tmpY to the plane conformed of X and Y
				myInfo.beforeX=tmpX-(tmpX*N)*N;
				myInfo.beforeY=CGAL::cross_product(myInfo.V,myInfo.beforeX);
				myInfo.beforeX=myInfo.beforeX/std::sqrt(myInfo.beforeX.squared_length());
				myInfo.beforeY=myInfo.beforeY/std::sqrt(myInfo.beforeY.squared_length());
			}

			for(double ss=0.5;ss<myInfo.DIV;ss++)
			{
				double s=((double)ss)/((double)myInfo.DIV);
				for(auto particle=myInfo.particles.begin();particle!=myInfo.particles.end();particle++)
				{
					Vector B=myInfo.beforeX*(*particle).x()+myInfo.beforeY*(*particle).y();
					int li, lj;
					Locate_type lt;

					//Project N
					double c1=-(B*myInfo.itrD->Normal)/(myInfo.V*myInfo.itrD->Normal);
					double c2=-(B*(myInfo.itrD+1)->Normal)/(myInfo.V*(myInfo.itrD+1)->Normal);
					Vector tmp1=B+c1*myInfo.V;
					Vector tmp2=B+c2*myInfo.V;
					Point D1=(*myInfo.itrC)+tmp1;
					Point D2=(*(myInfo.itrC+1))+tmp2;
					Point D(D2.x()*s+D1.x()*(1-s),D2.y()*s+D1.y()*(1-s),D2.z()*s+D1.z()*(1-s));
					if(myInfo.before==NULL)
						c = T.locate(D, lt,li,lj);
					else
						c=T.locate(D,lt,li,lj,myInfo.before);
					
					
					if(lt==Locate_type::CELL)
					{
							std::map<Delaunay::Cell_handle,__int64>::iterator it_c=index.find(c);
							__int64 d=it_c->second;
							cells[d]++;
							myInfo.before=c;
					}
					
					++numInterior;
				}
			}
	
	}, size);

	for(vector<Rad_branch>::iterator itr=data.begin();itr!=data.end();itr++)
	{
		Rad_branch a=*itr;
		const branch* _branch=boost::get<1>(a);
		delete(_branch);
	}
	data.clear();
	for(vector<eclipses*>::iterator itr=eclipseTree.begin();itr!=eclipseTree.end();itr++)
	{
		delete(*itr);
	}
	eclipseTree.clear();
	std::cout<<"start cell locate"<<endl;	
	__int64 N=0;

	std::cout<<endl;
	
	std::cout<<"T.number_of_finite_cells:"<<T.number_of_finite_cells()<<endl;
	std::cout<<"T.number_of_finite_facets:"<<T.number_of_finite_facets()<<endl;
	bool* bool_list=new bool[T.number_of_finite_cells()];

	N=0;
	double D=0.4;
	__int64 count=0;
	for(Delaunay::Finite_cells_iterator itr=T.finite_cells_begin();itr!=T.finite_cells_end();itr++,N++)
	{
		bool_list[N]=true;
		if(cells[N]<1)bool_list[N]=false;

	}
	std::cout << "releasing cells" << endl;
	delete [] cells;

	std::cout<<"start refine"<<endl;
	//int everything=1;
	//while(everything>0){
	//everything=0;
	std::cout<<"erase bubbles"<<endl;
	std::map<Cell_handle,__int64> _cells; 
	std::list<Cell_handle> cells1; 
	std::list<Cell_handle> cells2; 
	std::vector<std::list<Cell_handle>> cell_group;
	__int64 totalCount=0;
	N=0;
	for(auto itr=T.finite_cells_begin();itr!=T.finite_cells_end();itr++)
	{
		if(bool_list[index.find(itr)->second]==false)
		{
			_cells.insert(std::make_pair(itr,N));
			totalCount++;
			N++;
		}
	}
	__int64 NN=totalCount/20;
	if(NN=0)NN=1;
	cells2.push_back(_cells.begin()->first);
	_cells.erase(_cells.begin()->first);
	N=0;
	std::cout<<"totalCount:"<<totalCount<<endl;
	while(!_cells.empty())
	{
		std::list<Cell_handle> cells3; 
		while(!cells2.empty())
		{
			if(!_cells.empty())
			{
				for(auto itr=cells2.begin();itr!=cells2.end();itr++)
				{
					for(int i=0;i<4;i++)
					{
						Cell_handle nei=(*itr)->neighbor(i);
						std::map<Cell_handle,__int64>::iterator it=_cells.find(nei);
						if(it!=_cells.end())
						{
							cells1.push_back(nei);
							_cells.erase(nei);
							N++;
							if(((__int64)N/NN)*NN==N)std::cout<<"*";
						}
					}
				}
			}
			for(auto itr=cells2.begin();itr!=cells2.end();itr++)
			{
				cells3.push_back(*itr);
			}
			cells2.clear();
			if(!cells1.empty()){
				for(auto itr=cells1.begin();itr!=cells1.end();itr++)
				{
					cells2.push_back(*itr);
				}
				cells1.clear();
			}
		}
		if(!cells2.empty())
		{
			for(auto itr=cells2.begin();itr!=cells2.end();itr++)
			{
				cells3.push_back(*itr);
			}
		}
		if(!cells3.empty())
			cell_group.push_back(cells3);
		if(!_cells.empty())
		{
			cells2.push_back(_cells.begin()->first);
			_cells.erase(_cells.begin()->first);
		}
	}
	std::cout<<endl;
	std::cout<<cell_group.size()<<"groups found"<<endl;
	count=0;
	parallel_for_each(cell_group.begin(), cell_group.end(), [&index,&bool_list,&count](std::list<Cell_handle> itr)
	//for(auto itr=cell_group.begin();itr!=cell_group.end();itr++)
	{
		bool flag=false;
		for(auto itr2=(itr).begin();itr2!=(itr).end();itr2++)
		{
			for(int i=0;i<4;i++)
			{
				if(index.find((*itr2)->neighbor(i))==index.end())
				{
					flag=true;
				}
			}
			if(flag)break;
		}
		if(!flag)
		{
			for(auto itr2=(itr).begin();itr2!=(itr).end();itr2++)
			{
				bool_list[index.find(*itr2)->second]=true;
				//everything++;
				count++;
			}
		}
	});
	std::cout<<count<<"cells recovered"<<endl;
	
	std::cout<<"erase irregular incident vertices"<<endl;
	while(true)
	{
		N=0;
		NN=T.number_of_vertices()/20;
		if(NN==0)NN=1;
		totalCount=0;
		auto itr = T.vertices_begin();
		for(auto itr=T.vertices_begin();itr!=T.vertices_end();itr++,N++)
		{
			std::list<Cell_handle> _cells; 
			std::list<Cell_handle> cells; 
			std::list<Cell_handle> cells2; 
			std::list<Cell_handle> cells3; 
			T.incident_cells(itr,std::back_inserter(_cells)); 
			if (_cells.empty()){
				continue;
			}
			for(auto _itr=_cells.begin();_itr!=_cells.end();_itr++)//_itr was itr, possibly cause an error.
			{
				map<Cell_handle,__int64>::iterator itr2=index.find(*_itr);
				if(itr2!=index.end())
				{
					__int64 N=itr2->second;
					if(bool_list[N])
					{
						cells.push_back(*_itr);
					}
				}
			}
			if (cells.empty()){
				continue;
			}
			if (cells.size() == 1){
				continue;
			}
			cells2.push_back(*cells.begin());
			cells.pop_front();
			while(true)
			{
				int count=0;
				for(auto _itr=cells.begin();_itr!=cells.end();_itr++)
				{
					for(int i=0;i<4;i++)
					{
						Cell_handle nei=(*_itr)->neighbor(i);
						if(std::find(cells2.begin(),cells2.end(),nei)!=cells2.end())
						{
							cells3.push_back(*_itr);
							count++;
							break;
						}
					}
				}
				for(auto itr=cells3.begin();itr!=cells3.end();itr++)
				{
					cells.erase(std::find(cells.begin(),cells.end(),*itr));
					cells2.push_back(*itr);
				}
				cells3.clear();
				if(count==0)break;
			}
			if(cells.size()!=0){
				std::list<Cell_handle> toErase;
				if(cells.size()<cells2.size()) toErase=cells; else toErase=cells2;
				for(auto itr=toErase.begin();itr!=toErase.end();itr++)
				{
					__int64 N=index.find(*itr)->second;
					bool_list[N]=false;
					//everything++;
					totalCount++;
				}
			}
			if(((int)N/NN)*NN==N)std::cout<<"*";
		}
		std::cout<<endl;
		std::cout<<totalCount<<"cells removed!"<<endl;
		if(totalCount==0)break;
	}
	std::cout<<"erase irregular incident edges"<<endl;
	while(true)
	{
		N=0;
		NN=T.number_of_edges()/20;
		if(NN==0)NN=1;
		totalCount=0;
		for(auto itr=T.finite_edges_begin();itr!=T.finite_edges_end();itr++,N++)
		{
			std::list<Cell_handle> _cells; 
			std::list<Cell_handle> cells; 
			std::list<Cell_handle> cells2; 
			std::list<Cell_handle> cells3; 
			Delaunay::Cell_circulator circle=T.incident_cells(*itr); 
			Delaunay::Cell_circulator begin=circle;
			std::map<Cell_handle,__int64>::iterator it=index.find(circle);
			if(it==index.end())continue;
			bool flag=bool_list[it->second];
			int counter=0;
			std::vector<Delaunay::Cell_circulator> start;
			while(true)
			{
				circle++;
				it=index.find(circle);
				if(it==index.end())
				{
					if(flag==true)
					{
						flag=false;
						counter++;
					}
				}else if(bool_list[it->second]!=flag){
					flag=bool_list[it->second];
					if(flag)
					{
						start.push_back(circle);
					}
					counter++;
				}
				if(circle==begin)break;
			}
			if(counter==4)
			{
				std::vector<__int64> len;
				for(int i=0;i<2;i++)
				{
					Delaunay::Cell_circulator s=start[i];
				    __int64 L=0;
					while(true)
					{
						circle++;
						L++;
						it=index.find(circle);
						if(it==index.end())
						{
							break;
						}else if(bool_list[it->second]==false){
							break;
						}
					}
					len.push_back(L);
				}
				Delaunay::Cell_circulator toErase;
				if(len[0]<len[1]){
					toErase=start[0];
				}else
				{
					toErase=start[1];
				}
				while(true)
				{
					it=index.find(toErase);
					if(it==index.end())
					{
						break;
					}else if(bool_list[it->second]==false){
						break;
					}
					bool_list[it->second]=false;
					toErase++;
					totalCount++;
				}
			}
		}
		std::cout<<endl;
		std::cout<<totalCount<<"cells removed!"<<endl;
		if(totalCount==0)break;
	}
	
	std::cout<<"push and pop"<<endl;
	__int64 TT=0;
	while(true)
	{
		int num=0;
		N=0;
		for(Delaunay::Finite_cells_iterator itr=T.finite_cells_begin();itr!=T.finite_cells_end();itr++,N++)
		{
			if(bool_list[N]==false)
			{
				int count=0;
				for(int i=0;i<4;i++)
				{
					Delaunay::Cell_handle _neighbor=itr->neighbor(i);
					std::map<Delaunay::Cell_handle,__int64>::iterator it_N=index.find(_neighbor);
					if(it_N!=index.end())
					{
						if(bool_list[it_N->second])count++;
					}
				}
				if(count>2){
					bool_list[N]=true;
					//everything++;
					num++;
				}
			}

		}

		TT++;
		std::cout<<"refine:"<<TT<<", corrected:"<<num<<endl;
		if(num==0)break;
	}
	std::cout << "look up isolated tets." << endl;
	N = 0;
	__int64 num = 0;
	for (Delaunay::Finite_cells_iterator itr = T.finite_cells_begin(); itr != T.finite_cells_end(); itr++, N++)
	{
		if (bool_list[N] == true)
		{
			int count = 0;
			for (int i = 0; i<4; i++)
			{
				Delaunay::Cell_handle _neighbor = itr->neighbor(i);
				std::map<Delaunay::Cell_handle, __int64>::iterator it_N = index.find(_neighbor);
				if (it_N != index.end())
				{
					if (!bool_list[it_N->second])count++;
				}
				else
				{
					count++;
				}
			}
			if (count==4){
				bool_list[N] = false;
				//everything++;
				num++;
			}
		}

	}

	std::cout << "found "<<num<<" isolated tets." << endl;

	NN=T.number_of_finite_facets()/20;
	if(NN<1)NN=1;
	N=0;
	std::vector<Delaunay::Facet> facet_list;
	std::cout<<"count up boundary facets"<<endl;
	for(Delaunay::Finite_cells_iterator itr=T.finite_cells_begin();itr!=T.finite_cells_end();itr++,N++)
	{
		if(bool_list[N])
		{
			for(int i=0;i<4;i++)
			{
				Delaunay::Cell_handle nei=itr->neighbor(i);
				
				std::map<Delaunay::Cell_handle,__int64>::iterator itr_c=index.find(nei);
				if(itr_c==index.end())
				{
					facet_list.push_back(Delaunay::Facet(itr,i));
				}else
				{
					__int64 N2=itr_c->second;
					if(!bool_list[N2])
						facet_list.push_back(Delaunay::Facet(itr,i));
				}
			}
		}
		if(((int)N/NN)*NN==N)std::cout<<"*";
	}
	
	std::cout<<endl;
	std::cout << "releasing bool_list" << endl;
	delete [] bool_list;
	std::cout << "T.number_of_finite_facets:" << T.number_of_finite_facets() << endl;
	std::cout<<"start writing file"<<"["<<filename1<<"]"<<endl;
	
	N=0;
	NN=T.number_of_finite_cells()/20;
	if(NN==0)NN=1;
	/*
	for(auto v :exterior)
	{
		ofs3<<v.x()<<" , "<<v.y()<<" , "<<v.z()<<endl;
	}*/
	std::cout<<endl;
	N=0;
	num=0;
	NN=facet_list.size()/20;
	if(NN==0)NN=1;
	
	std::map<Delaunay::Vertex_handle,__int64> vIndex;

	for(auto itr=facet_list.begin();itr!=facet_list.end();itr++,num++)
	{
		for(int i=0;i<4;i++)
		{
			if(i!=itr->second)
			{
				Delaunay::Vertex_handle handle=itr->first->vertex(i);				
				std::map<Delaunay::Vertex_handle,__int64>::iterator pair=vIndex.find(handle);
				if(pair==vIndex.end())
				{
					Delaunay::Point P=handle->point();
					double x=P.x(),y=P.y(),z=P.z();
					ofs<<x<<" , "<<y<<" , "<<z<<endl;
					vIndex.insert(std::make_pair(handle,N));
					ofs2<<N<<" ";
					N++;
				}else
				{
					ofs2<<pair->second<<" ";
				}
			}
		}
		ofs2<<endl;
		if(((int)num/NN)*NN==num)std::cout<<"*";
	}
	
	std::cout<<endl;
	ofs.close();
	ofs2.close();
	//ofs3.close();
	
	//ofs4.close();
	
	
	std::cout << "Press Return To Exit...";
	std::cin.get();




  return 0;
}
