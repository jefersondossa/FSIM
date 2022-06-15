#pragma once
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
















Surface::Surface() {}

Surface::Surface(const int& index, const std::string& name, LineLoop* lineLoop)
{
	index_ = index;
	name_ = name;
	lineLoop_ = lineLoop;
}

Surface::~Surface() {}

int Surface::getIndex()
{
	return index_;
}

std::string Surface::getName()
{
	return name_;
}

LineLoop* Surface::getLineLoop()
{
	return lineLoop_;
}

Material* Surface::getMaterial()
{
	return material_;
}

// std::vector<ParametricSurfaceElement*> Surface::getParametricElements()
// {
// 	return parametricElements_;
// }

std::vector<ElementMesh*> Surface::getElements()
{
	return elements_;
}

std::vector<NodeMesh*> Surface::getNodes()
{
	return nodes_;
}

void Surface::setMaterial(Material* material)
{
	material_ = material;
}

void Surface::addElement(ElementMesh* element)
{
	elements_.push_back(element);
}

// void Surface::addParametricElement(ParametricSurfaceElement* element)
// {
// 	parametricElements_.push_back(element);
// }

std::string Surface::getGmshCode()
{
	std::stringstream text;
	text << name_ << " = news; Surface(" << name_ << ") = {" << lineLoop_->getName() << "}; Physical Surface('" << name_ << "') = {" << name_ << "};\n//\n";
	return text.str();
}

void Surface::addNodes(const std::vector<NodeMesh*>& nodes)
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

void Surface::clearAllElements()
{
	elements_.clear();
	// parametricElements_.clear();
	std::vector<ElementMesh*>().swap(elements_);
	// std::vector<ParametricSurfaceElement*>().swap(parametricElements_);
}

void Surface::clearAllNodes()
{
	nodes_.clear();
	std::vector<NodeMesh*>().swap(nodes_);
}

void Surface::removeElement(const int& index)
{
	elements_.erase(elements_.begin() + index);
	// parametricElements_.erase(parametricElements_.begin() + index);
}
