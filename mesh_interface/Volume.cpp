#include "Volume.h"





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

std::vector<NodeMesh*> Volume::NodeVec()
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
