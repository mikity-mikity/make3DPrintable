#ifndef _HALFEDGE_H_
#define _HALFEDGE_H_
#include <vector>
#include "face.h"
#include "vertex.h"
using namespace std;
namespace GeometryProcessing
{
	class face;
	class vertex;
    class halfedge
    {
	public:
		vertex *P;
        face *owner;
        halfedge *pair, *next, *prev;
        bool isNaked();
		bool isBoundary();
        halfedge(vertex *_P);
    };
}
#endif