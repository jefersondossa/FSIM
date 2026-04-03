#ifndef VTUREADER_H
#define VTUREADER_H

#include <string>
#include <vector>

namespace VTUReader
{
using Grid = std::vector<std::vector<double>>;

Grid read_field_from_vtu(const std::string& filename, const std::string& field_name);
}

#endif
