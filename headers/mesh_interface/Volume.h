#pragma once

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






Volume::Volume() {}

Volume::Volume(const int& index, const std::string& name, SurfaceLoop* surfaceLoop)
{
	index_ = index;
	name_ = name;
	surfaceLoop_ = surfaceLoop;
}

Volume::~Volume() {}

int Volume::getIndex()
{
	return index_;
}

std::string Volume::getName()
{
	return name_;
}

SurfaceLoop* Volume::getSurfaceLoop()
{
	return surfaceLoop_;
}

Material* Volume::getMaterial()
{
	return material_;
}

// std::vector<ParametricVolumeElement*> Volume::getParametricElements()
// {
// 	return parametricElements_;
// }

std::vector<ElementMesh*> Volume::getElements()
{
	return elements_;
}

std::vector<NodeMesh*> Volume::getNodes()
{
	return nodes_;
}

void Volume::setMaterial(Material* material)
{
	material_ = material;
}

void Volume::addElement(ElementMesh* element)
{
	elements_.push_back(element);
}

// void Volume::addParametricElement(ParametricVolumeElement* element)
// {
// 	parametricElements_.push_back(element);
// }

void Volume::addNodes(const std::vector<NodeMesh*>& nodes)
{
	for (NodeMesh* node1 : nodes)
	{
		bool notDuplicate = true;
		for (NodeMesh* node2 : nodes_)
		{
			if (node1->getIndex() == node2->getIndex())
			{
				notDuplicate = false;
				break;
			}
		}
		if (notDuplicate)
			nodes_.push_back(node1);
	}
}

std::string Volume::getGmshCode()
{
	std::stringstream text;
	text << name_ << " = newv; Volume(" << name_ << ") = {" << surfaceLoop_->getName() << "}; Physical Volume('" << name_ << "') = {" << name_ << "};\n//\n";
	return text.str();
}
