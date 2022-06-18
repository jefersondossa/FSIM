#ifndef SURFACE_H
#define SURFACE_H

#include "LineLoop.h"
#include "Material.h"
// #include "ParametricSurfaceElement.h"
#include <algorithm>

class Surface
{
public:
	Surface();

	Surface(const int& index, const std::string& name, LineLoop* lineLoop);

	~Surface();

	int getIndex();

	std::string getName();

	LineLoop* getLineLoop();

	Material* getMaterial();

	// std::vector<ParametricSurfaceElement*> getParametricElements();

	std::vector<ElementMesh*> getElements();

	std::vector<NodeMesh*> getNodes();

	void setMaterial(Material* material);

	void addElement(ElementMesh* element);

	// void addParametricElement(ParametricSurfaceElement* element);

	virtual std::string getGmshCode();

	void addNodes(const std::vector<NodeMesh*>& nodes);

	void clearAllElements();

	void clearAllNodes();

	void removeElement(const int& index);

protected:
	int index_;
	std::string name_;
	LineLoop* lineLoop_;
	Material* material_;
	// std::vector<ParametricSurfaceElement*> parametricElements_;
	std::vector<ElementMesh*> elements_;
	std::vector<NodeMesh*> nodes_;
};

#endif










