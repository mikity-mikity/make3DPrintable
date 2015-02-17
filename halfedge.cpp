#include "halfedge.h"
#include "vertex.h"
namespace GeometryProcessing{

	bool halfedge::isNaked()
	{
		if(pair==NULL)return true;
		return false;
	}
	halfedge::halfedge(vertex *_P)
	{
		P = _P;
		if (_P->hf_begin == NULL) _P->hf_begin = this;
		owner=NULL;
		pair=NULL;
		next=NULL;
		prev=NULL;
	}
}
