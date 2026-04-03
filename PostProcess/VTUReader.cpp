#include "VTUReader.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace
{
struct DataArray
{
    std::string name;
    int number_of_components = 1;
    std::string values;
};

struct Point
{
    double x = 0.0;
    double y = 0.0;
};

std::string read_text_file(const std::string& filename)
{
    std::ifstream input(filename);
    if (!input) {
        throw std::runtime_error("failed to open file: " + filename);
    }

    std::ostringstream buffer;
    buffer << input.rdbuf();
    return buffer.str();
}

std::unordered_map<std::string, std::string> parse_attributes(const std::string& tag)
{
    std::unordered_map<std::string, std::string> attributes;
    std::size_t pos = tag.find_first_of(" \t\r\n");
    if (pos == std::string::npos) {
        return attributes;
    }

    while (true) {
        pos = tag.find_first_not_of(" \t\r\n", pos);
        if (pos == std::string::npos || tag[pos] == '>' || tag[pos] == '/') {
            break;
        }

        const auto eq = tag.find('=', pos);
        if (eq == std::string::npos) {
            break;
        }

        const auto key_end = tag.find_last_not_of(" \t\r\n", eq - 1);
        if (key_end == std::string::npos || key_end < pos) {
            break;
        }

        const auto quote_begin = tag.find('"', eq + 1);
        if (quote_begin == std::string::npos) {
            break;
        }

        const auto quote_end = tag.find('"', quote_begin + 1);
        if (quote_end == std::string::npos) {
            break;
        }

        attributes[tag.substr(pos, key_end - pos + 1)] = tag.substr(quote_begin + 1, quote_end - quote_begin - 1);
        pos = quote_end + 1;
    }

    return attributes;
}

std::string extract_tag_block(const std::string& xml, const std::string& tag_name)
{
    const auto open = xml.find("<" + tag_name);
    if (open == std::string::npos) {
        throw std::runtime_error("missing <" + tag_name + "> section");
    }

    const auto open_end = xml.find('>', open);
    if (open_end == std::string::npos) {
        throw std::runtime_error("malformed <" + tag_name + "> tag");
    }

    const auto close = xml.find("</" + tag_name + ">", open_end + 1);
    if (close == std::string::npos) {
        throw std::runtime_error("missing </" + tag_name + "> section");
    }

    return xml.substr(open_end + 1, close - open_end - 1);
}

std::vector<DataArray> extract_data_arrays(const std::string& block)
{
    std::vector<DataArray> arrays;
    std::size_t pos = 0;

    while (true) {
        const auto open = block.find("<DataArray", pos);
        if (open == std::string::npos) {
            break;
        }

        const auto open_end = block.find('>', open);
        if (open_end == std::string::npos) {
            throw std::runtime_error("malformed <DataArray> tag");
        }

        const auto close = block.find("</DataArray>", open_end + 1);
        if (close == std::string::npos) {
            throw std::runtime_error("missing </DataArray> tag");
        }

        const auto tag = block.substr(open, open_end - open + 1);
        const auto attributes = parse_attributes(tag);

        DataArray array;
        if (const auto it = attributes.find("Name"); it != attributes.end()) {
            array.name = it->second;
        }
        if (const auto it = attributes.find("NumberOfComponents"); it != attributes.end()) {
            array.number_of_components = std::stoi(it->second);
        }
        array.values = block.substr(open_end + 1, close - open_end - 1);
        arrays.push_back(std::move(array));

        pos = close + std::string("</DataArray>").size();
    }

    return arrays;
}

std::vector<double> parse_doubles(const std::string& text)
{
    std::vector<double> values;
    std::istringstream input(text);
    double value = 0.0;
    while (input >> value) {
        values.push_back(value);
    }
    return values;
}

std::string normalize_field_name(const std::string& name)
{
    std::string normalized;
    normalized.reserve(name.size());

    for (unsigned char ch : name) {
        if (std::isalnum(ch)) {
            normalized.push_back(static_cast<char>(std::tolower(ch)));
        }
    }

    return normalized;
}

double coordinate_tolerance(const std::vector<double>& coordinates)
{
    if (coordinates.empty()) {
        return 1e-9;
    }

    auto [min_it, max_it] = std::minmax_element(coordinates.begin(), coordinates.end());
    const double scale = std::max({1.0, std::abs(*min_it), std::abs(*max_it), std::abs(*max_it - *min_it)});
    return 1e-9 * scale;
}

std::vector<double> unique_sorted_coordinates(std::vector<double> coordinates)
{
    std::sort(coordinates.begin(), coordinates.end());

    const double tol = coordinate_tolerance(coordinates);
    std::vector<double> unique_values;

    for (double coordinate : coordinates) {
        if (unique_values.empty() || std::abs(coordinate - unique_values.back()) > tol) {
            unique_values.push_back(coordinate);
        }
    }

    return unique_values;
}

std::size_t find_coordinate_index(const std::vector<double>& unique_values, double coordinate)
{
    const double tol = coordinate_tolerance(unique_values);
    auto it = std::lower_bound(unique_values.begin(), unique_values.end(), coordinate);

    if (it != unique_values.end() && std::abs(*it - coordinate) <= tol) {
        return static_cast<std::size_t>(std::distance(unique_values.begin(), it));
    }

    if (it != unique_values.begin()) {
        --it;
        if (std::abs(*it - coordinate) <= tol) {
            return static_cast<std::size_t>(std::distance(unique_values.begin(), it));
        }
    }

    throw std::runtime_error("point coordinate does not match the reconstructed grid");
}
}

namespace VTUReader
{
Grid read_field_from_vtu(const std::string& filename, const std::string& field_name)
{
    try {
        const auto xml = read_text_file(filename);

        const auto points_block = extract_tag_block(xml, "Points");
        const auto point_arrays = extract_data_arrays(points_block);
        if (point_arrays.empty()) {
            throw std::runtime_error("no point coordinates found");
        }

        const auto point_values = parse_doubles(point_arrays.front().values);
        const int point_components = point_arrays.front().number_of_components;
        if (point_components < 2) {
            throw std::runtime_error("point coordinates must have at least two components");
        }
        if (point_values.size() % static_cast<std::size_t>(point_components) != 0) {
            throw std::runtime_error("point coordinate array has an invalid size");
        }

        std::vector<Point> points;
        points.reserve(point_values.size() / static_cast<std::size_t>(point_components));
        for (std::size_t i = 0; i < point_values.size(); i += static_cast<std::size_t>(point_components)) {
            points.push_back({point_values[i], point_values[i + 1]});
        }

        const auto point_data_block = extract_tag_block(xml, "PointData");
        const auto point_data_arrays = extract_data_arrays(point_data_block);

        const DataArray* field_array = nullptr;
        for (const auto& array : point_data_arrays) {
            if (array.name == field_name) {
                field_array = &array;
                break;
            }
        }
        if (field_array == nullptr) {
            const auto normalized_name = normalize_field_name(field_name);
            for (const auto& array : point_data_arrays) {
                if (normalize_field_name(array.name) == normalized_name) {
                    field_array = &array;
                    break;
                }
            }
        }
        if (field_array == nullptr) {
            throw std::runtime_error("field not found in PointData: " + field_name);
        }
        if (field_array->number_of_components != 1) {
            throw std::runtime_error("field is not scalar: " + field_array->name);
        }

        const auto field_values = parse_doubles(field_array->values);
        if (field_values.size() != points.size()) {
            throw std::runtime_error("field size does not match the number of points");
        }

        std::vector<double> xs;
        std::vector<double> ys;
        xs.reserve(points.size());
        ys.reserve(points.size());
        for (const auto& point : points) {
            xs.push_back(point.x);
            ys.push_back(point.y);
        }

        const auto unique_x = unique_sorted_coordinates(std::move(xs));
        const auto unique_y = unique_sorted_coordinates(std::move(ys));
        if (unique_x.size() * unique_y.size() != points.size()) {
            throw std::runtime_error("points do not form a rectangular grid");
        }

        // Row-major layout: grid[row][col] maps to grid[y][x].
        Grid grid(unique_y.size(), std::vector<double>(unique_x.size(), 0.0));
        std::vector<std::vector<bool>> filled(unique_y.size(), std::vector<bool>(unique_x.size(), false));

        for (std::size_t point_index = 0; point_index < points.size(); ++point_index) {
            const std::size_t row = find_coordinate_index(unique_y, points[point_index].y);
            const std::size_t col = find_coordinate_index(unique_x, points[point_index].x);

            if (filled[row][col]) {
                throw std::runtime_error("duplicate point found at the same grid position");
            }

            grid[row][col] = field_values[point_index];
            filled[row][col] = true;
        }

        for (std::size_t row = 0; row < filled.size(); ++row) {
            for (std::size_t col = 0; col < filled[row].size(); ++col) {
                if (!filled[row][col]) {
                    throw std::runtime_error("grid reconstruction left missing points");
                }
            }
        }

        return grid;
    } catch (const std::exception& error) {
        std::cerr << "read_field_from_vtu error: " << error.what() << '\n';
        return {};
    }
}
}
