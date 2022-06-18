#ifndef VOLUME_H
#define VOLUME_H

#include "SurfaceLoop.h"
#include "Material.h"
// #include "ParametricVolumeElement.h"

class Volume
{
public:
	Volume();

	Volume(const int& index, const std::string& name, SurfaceLoop* surfaceLoop);

	~Volume();

	int getIndex();

	std::string getName();

	SurfaceLoop* getSurfaceLoop();

	Material* getMaterial();

	// std::vector<ParametricVolumeElement*> getParametricElements();

	std::vector<ElementMesh*> getElements();

	std::vector<NodeMesh*> getNodes();

	void setMaterial(Material* material);

	void addElement(ElementMesh* element);

	// void addParametricElement(ParametricVolumeElement* element);

	void addNodes(const std::vector<NodeMesh*>& nodes);

	std::string getGmshCode();

private:
	int index_;
	std::string name_;
	SurfaceLoop* surfaceLoop_;
	Material* material_;
	std::vector<ElementMesh*> elements_;
	std::vector<NodeMesh*> nodes_;
	// std::vector<ParametricVolumeElement*> parametricElements_;
};


#endif