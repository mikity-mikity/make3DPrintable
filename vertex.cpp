#include"vertex.h"
#include"halfedge.h"
namespace GeometryProcessing
{
	bool vertex::isNaked()
	{
		if(hf_begin == hf_end) return false;
		return true;
	}
	vertex::vertex(int _N)
	{
		N = _N;
		hf_begin=NULL;
		hf_end=NULL;
	}
	bool vertex::isInner()
	{
		return onering[0] == onering[onering.size() - 1]->next;
	}
	bool vertex::isBoundary()
	{
		return onering[0] != onering[onering.size() - 1]->next;
	}
}