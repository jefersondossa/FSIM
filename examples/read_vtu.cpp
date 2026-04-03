#include <VTUReader.h>

int main()
{
    const auto filename = "cantilever_2d_beam0.vtu";
    const auto field_name = "weight_function";
    auto data = VTUReader::read_field_from_vtu(filename, field_name);
    return 0;
}
