#include <vector>
#include<Eigen/Dense>
#include<Eigen/Sparse>

using namespace std;
using namespace Eigen;
namespace GeometryProcessing
{
	struct MeshFace
	{
	public:
		int A,B,C;
	};
	struct Point3d{
	public:
		double X,Y,Z;
	};
	class Mesh
	{
	public:
		vector<MeshFace> Faces;
		vector<Point3d> Vertices;
	};
	class vertex
    {
	public:
		int N;
        vector<halfedge*> star;
        vector<halfedge*> onering;
        halfedge* hf_begin;
        halfedge* hf_end;
        bool isNaked()
        {
            if(hf_begin == hf_end) return false;
			return true;
        }
        vertex(int _N)
        {
            N = _N;
        }
        bool isInner()
        {
            return onering[0] == onering[onering.size() - 1]->next;
        }
        bool isBoundary()
        {
            return onering[0] != onering[onering.size() - 1]->next;
        }
    };

    class halfedge
    {
	public:
		vertex *P;
        face *owner;
        halfedge *pair, *next, *prev;
        bool isNaked()
        {
            if(pair)return true;
			return false;
        }
        halfedge(vertex *_P)
        {
            P = _P;
            if (_P->hf_begin == NULL) _P->hf_begin = this;
        }
    };
    class face
    {
	public:
        int N;
        int corner[3];
        halfedge *firsthalfedge;
		face(int _N, int A,int B,int C)
        {
            corner[0]=A;
			corner[1]=B;
			corner[2]=C;
            N = _N;
        }
    };

    class MeshStructure
    {
	private:
		enum orient
        {
            unknown, positive, negative
        };

        //to get boundary chain
        //boundaryStart->hf->next->P->hf->next->P->....
	public:
		vertex *boundaryStart;
        vector<vertex*> vertices;
        vector<face*> faces;
        vector<halfedge*> halfedges;
        vector<vertex*> innerVertices;
        vector<vertex*> outerVertices;  
	private:
		//halfedge** __halfedgeTable;
        Eigen::SparseMatrix<halfedge*> __halfedgeTable;
		Eigen::SparseMatrix<vector<face*>*> _faceTable;
        orient* __orientation;
	public:
		vector<halfedge*> edges()
        {
            auto res = vector<halfedge*>();
            for (auto e : halfedges)
            {
                if (e->isNaked())
                {
                    res.push_back(e);
                }
                else
                {
                    if (e->P->N < e->next->P->N)
                        res.push_back(e);
                }
            }
            return res;
        }
        int nVertices()
        {
            return vertices.size();
        }
        int nFaces()
        {
            return faces.size();
        }
	private:
		MeshStructure()
        {
            this->Clear();
        }
        void Construct(Mesh val)
        {
            int _nVertices = val.Vertices.size();
            int _nFaces = val.Faces.size();

            __orientation = new orient[_nFaces];
            _faceTable = SparseMatrix<vector<face*>*>(_nVertices, _nVertices);
            __halfedgeTable = SparseMatrix<halfedge*>(_nVertices, _nVertices);

            for (int i = 0; i < _nFaces; i++)
            {
                __orientation[i] = orient::unknown;
            }

            for (int i = 0; i < _nVertices; i++)
            {
                auto _v = new vertex(i);
                vertices.push_back(_v);
            }

            for (int i = 0; i < _nFaces; i++)
            {
                auto f = val.Faces[i];
                auto _f = new face(i, f.A, f.B, f.C);
                faces.push_back(_f);
                faceTableAdd(_f);
            }
            //Recursive
            halfEdgeAdd(faces[0]);
            //find pairs
            for (auto h : halfedges)
            {
                int i = h->P->N;
                int j = h->next->P->N;
				//if (__halfedgeTable.coeff(i, j) != NULL) throw new ArgumentOutOfRangeException(";)");
                __halfedgeTable.coeffRef(i, j) = h;
            }
            for (auto h : halfedges)
            {
                int i = h->P->N;
                int j = h->next->P->N;
                //if boundary edge...
                if (__halfedgeTable.coeff(j, i) == NULL)
                {
                    h->pair = NULL;
                }
                else
                {
                    h->pair = __halfedgeTable.coeffRef(j, i);
                }
            }
            //post process to find boundary vertices
            for (auto v : vertices)
            {
                auto h = v->hf_begin;
                v->hf_end = h;
                do
                {
                    if (h->prev->isNaked())
                    {
                        v->hf_end = h->prev;
                        while (!h->isNaked())
                        {
                            h = h->pair->next;
                        }
                        v->hf_begin = h;
                        if (this->boundaryStart ==NULL) this->boundaryStart = v;
                        break;
                    }
                    h = h->prev->pair;
                } while (h != v->hf_begin);
            }

            //post process to create stars
            for (auto v : vertices)
            {
                auto h = v->hf_begin;
                v->star.clear();
                do
                {
                    v->star.push_back(h);
                    if (h->prev->isNaked()) break;
                    h = h->prev->pair;
                } while (h != v->hf_begin);
            }
            //post process to create onering
            for (auto v : vertices)
            {
                auto h = v->hf_begin;
                v->onering.clear();
                do
                {
                    do
                    {
                        h = h->next;
                        v->onering.push_back(h);
                    } while (h->next->next->P != v);
                    if (h->next->isNaked()) break;
                    h = h->next->pair;
                } while (h != v->hf_begin);
            }
            //post process to split the vertices into inner and outer.
            innerVertices.clear();
            outerVertices.clear();
            for (auto v : vertices)
            {
                if (v->hf_begin->isNaked()) outerVertices.push_back(v); else innerVertices.push_back(v);
            }
        }
		private:
			void halfEdgeAdd(face *f)
			{
				auto _o = orient::unknown;
				for (int i = 0; i < 3; i++)
				{
					int I = f->corner[i];
					int J = 0;
					if(i == 2) J=f->corner[0];else J=f->corner[i + 1];
					if (_faceTable.coeffRef(I, J)->size() == 2)
					{
						if (_faceTable.coeffRef(I, J)->at(0) == f)
						{
							if (__orientation[_faceTable.coeffRef(I, J)->at(1)->N] != orient::unknown)
							{
								_o = orient::unknown;
								if(__orientation[_faceTable.coeffRef(I, J)->at(1)->N] == orient::positive)
									_o=orient::negative;else orient::positive;
							}
						}
						if (_faceTable.coeffRef(I, J)->at(1) == f)
						{
							if (__orientation[_faceTable.coeffRef(I, J)->at(0)->N] != orient::unknown)
							{
								_o = orient::unknown;
								if(__orientation[_faceTable.coeffRef(I, J)->at(0)->N] == orient::positive)
									_o=orient::negative;else _o=orient::positive;
							}
						}
					}
					else
					{
						if (_faceTable.coeff(J, I) != NULL)
						{
							if (__orientation[_faceTable.coeffRef(J, I)->at(0)->N] != orient::unknown)
							{
								_o = __orientation[_faceTable.coeffRef(J, I)->at(0)->N];
							}
						}
					}
				}
				__orientation[f->N] = orient::unknown;
				if(_o == orient::unknown)
					__orientation[f->N]=orient::positive;else __orientation[f->N]=_o;
				//register a halfedge
				if (__orientation[f->N] == orient::positive)
				{
					auto he1 = new halfedge(vertices[f->corner[0]]);
					auto he2 = new halfedge(vertices[f->corner[1]]);
					auto he3 = new halfedge(vertices[f->corner[2]]);
					halfedges.push_back(he1);
					halfedges.push_back(he2);
					halfedges.push_back(he3);
					he1->prev = he3; he1->next = he2; he1->owner = f;
					he2->prev = he1; he2->next = he3; he2->owner = f;
					he3->prev = he2; he3->next = he1; he3->owner = f;
					f->firsthalfedge = he1;
				}

				if (__orientation[f->N] == orient::negative)
				{
					auto he1 = new halfedge(vertices[f->corner[2]]);
					auto he2 = new halfedge(vertices[f->corner[1]]);
					auto he3 = new halfedge(vertices[f->corner[0]]);
					halfedges.push_back(he1);
					halfedges.push_back(he2);
					halfedges.push_back(he3);
					he1->prev = he3; he1->next = he2; he1->owner = f;
					he2->prev = he1; he2->next = he3; he2->owner = f;
					he3->prev = he2; he3->next = he1; he3->owner = f;
					f->firsthalfedge = he1;
				}


				//list up neighbors that are not oriented
				for (int i = 0; i < 3; i++)
				{
					int I = f->corner[i];
					int J = 0;
					if(i == 2) J=f->corner[0];else J=f->corner[i + 1];
					if (_faceTable.coeffRef(I, J)->size() == 2)
					{
						if (_faceTable.coeffRef(I, J)->at(0) == f)
						{
							if (__orientation[_faceTable.coeffRef(I, J)->at(1)->N] == orient::unknown)
							{
								halfEdgeAdd(_faceTable.coeffRef(I, J)->at(1));
							}
						}
						if (_faceTable.coeffRef(I, J)->at(1) == f)
						{
							if (__orientation[_faceTable.coeffRef(I, J)->at(0)->N] == orient::unknown)
							{
								halfEdgeAdd(_faceTable.coeffRef(I, J)->at(0));
							}
						}
					}
					else
					{
						if (_faceTable.coeff(J, I) != NULL)
						{
							if (__orientation[_faceTable.coeffRef(J, I)->at(0)->N] == orient::unknown)
							{
								halfEdgeAdd(_faceTable.coeffRef(J, I)->at(0));
							}
						}
					}
				}
			}
		private:
			void faceTableAdd(int i, int j, face* f)
			{
				if (_faceTable.coeff(i, j)==NULL)
				{
					_faceTable.coeffRef(i, j) = new vector<face*>();
				}
				_faceTable.coeffRef(i, j)->push_back(f);
			}
			void faceTableAdd(face* f)
			{
				for (int i = 0; i < 3; i++)
				{
					int I = f->corner[i];
					int J=0;
					if(i == 2) J=f->corner[0];else J= f->corner[i + 1];
					faceTableAdd(I, J, f);
				}
			}
		public:
			static MeshStructure* CreateFrom(Mesh val)
			{
				MeshStructure* ret;
				ret->Construct(val);
				return ret;
			}
		public:
			void Clear()
			{
				vertices.clear();
				faces.clear();
				halfedges.clear();
				innerVertices.clear();
				outerVertices.clear();
				boundaryStart = NULL;
			}
    };
}