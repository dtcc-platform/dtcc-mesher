#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "dtcc_mesher/dtcc_mesher.h"
#include "dtcc_mesher/dtcc_mesher_version.h"

namespace py = pybind11;

static std::vector<tm_point> tm_points_from_array(const py::array_t<double, py::array::c_style | py::array::forcecast> &array)
{
    auto view = array.unchecked<2>();
    std::vector<tm_point> points(static_cast<size_t>(view.shape(0)));

    for (py::ssize_t i = 0; i < view.shape(0); ++i) {
        points[static_cast<size_t>(i)].x = view(i, 0);
        points[static_cast<size_t>(i)].y = view(i, 1);
    }

    return points;
}

static std::vector<tm_segment> tm_segments_from_array(const py::array_t<std::uint32_t, py::array::c_style | py::array::forcecast> &array)
{
    auto view = array.unchecked<2>();
    std::vector<tm_segment> segments(static_cast<size_t>(view.shape(0)));

    for (py::ssize_t i = 0; i < view.shape(0); ++i) {
        segments[static_cast<size_t>(i)].a = view(i, 0);
        segments[static_cast<size_t>(i)].b = view(i, 1);
    }

    return segments;
}

static py::array_t<double> tm_points_to_numpy(const tm_mesh &mesh)
{
    py::array_t<double> array({static_cast<py::ssize_t>(mesh.num_points), static_cast<py::ssize_t>(2)});
    auto view = array.mutable_unchecked<2>();

    for (size_t i = 0; i < mesh.num_points; ++i) {
        view(static_cast<py::ssize_t>(i), 0) = mesh.points[i].x;
        view(static_cast<py::ssize_t>(i), 1) = mesh.points[i].y;
    }

    return array;
}

static py::array_t<std::uint32_t> tm_triangles_to_numpy(const tm_mesh &mesh)
{
    py::array_t<std::uint32_t> array({static_cast<py::ssize_t>(mesh.num_triangles), static_cast<py::ssize_t>(3)});
    auto view = array.mutable_unchecked<2>();

    for (size_t i = 0; i < mesh.num_triangles; ++i) {
        view(static_cast<py::ssize_t>(i), 0) = mesh.triangles[i * 3 + 0];
        view(static_cast<py::ssize_t>(i), 1) = mesh.triangles[i * 3 + 1];
        view(static_cast<py::ssize_t>(i), 2) = mesh.triangles[i * 3 + 2];
    }

    return array;
}

static py::array_t<std::uint32_t> tm_segments_to_numpy(const tm_mesh &mesh)
{
    py::array_t<std::uint32_t> array({static_cast<py::ssize_t>(mesh.num_segments), static_cast<py::ssize_t>(2)});
    auto view = array.mutable_unchecked<2>();

    for (size_t i = 0; i < mesh.num_segments; ++i) {
        view(static_cast<py::ssize_t>(i), 0) = mesh.segments[i].a;
        view(static_cast<py::ssize_t>(i), 1) = mesh.segments[i].b;
    }

    return array;
}

static py::dict tm_summary_to_dict(const tm_quality_summary &summary)
{
    py::dict result;

    result["point_count"] = py::int_(summary.point_count);
    result["input_point_count"] = py::int_(summary.input_point_count);
    result["steiner_point_count"] = py::int_(summary.steiner_point_count);
    result["segment_split_point_count"] = py::int_(summary.segment_split_point_count);
    result["triangle_split_point_count"] = py::int_(summary.triangle_split_point_count);
    result["protected_corner_count"] = py::int_(summary.protected_corner_count);
    result["exempt_triangle_count"] = py::int_(summary.exempt_triangle_count);
    result["triangle_count"] = py::int_(summary.triangle_count);
    result["area_min"] = py::float_(summary.area_min);
    result["area_mean"] = py::float_(summary.area_mean);
    result["area_max"] = py::float_(summary.area_max);
    result["min_angle_deg_min"] = py::float_(summary.min_angle_deg_min);
    result["min_angle_deg_mean"] = py::float_(summary.min_angle_deg_mean);
    result["min_angle_deg_max"] = py::float_(summary.min_angle_deg_max);
    result["edge_ratio_min"] = py::float_(summary.edge_ratio_min);
    result["edge_ratio_mean"] = py::float_(summary.edge_ratio_mean);
    result["edge_ratio_max"] = py::float_(summary.edge_ratio_max);
    result["radius_edge_ratio_min"] = py::float_(summary.radius_edge_ratio_min);
    result["radius_edge_ratio_mean"] = py::float_(summary.radius_edge_ratio_mean);
    result["radius_edge_ratio_max"] = py::float_(summary.radius_edge_ratio_max);
    result["count_min_angle_lt_20"] = py::int_(summary.count_min_angle_lt_20);
    result["count_min_angle_lt_30"] = py::int_(summary.count_min_angle_lt_30);
    return result;
}

static py::dict tm_generate_raw(
    const py::array_t<double, py::array::c_style | py::array::forcecast> &points_array,
    py::object segments_object,
    py::object holes_object,
    double min_angle_deg,
    bool enable_refinement,
    bool use_offcenters,
    bool verbose,
    int acute_mode,
    py::object protect_angle_object,
    std::size_t max_refinement_steps,
    std::size_t max_protection_levels
)
{
    std::vector<tm_point> points;
    std::vector<tm_segment> segments;
    std::vector<tm_point> holes;
    tm_domain domain{};
    tm_options options;
    tm_mesh mesh{};
    tm_quality_summary summary{};
    tm_error error{};
    tm_status status;

    if (points_array.ndim() != 2 || points_array.shape(1) != 2) {
        throw std::invalid_argument("points must have shape (N, 2)");
    }

    points = tm_points_from_array(points_array);
    domain.points = points.empty() ? nullptr : points.data();
    domain.num_points = points.size();

    if (!segments_object.is_none()) {
        auto segments_array = py::cast<py::array_t<std::uint32_t, py::array::c_style | py::array::forcecast>>(segments_object);
        if (segments_array.ndim() != 2 || segments_array.shape(1) != 2) {
            throw std::invalid_argument("segments must have shape (M, 2)");
        }
        segments = tm_segments_from_array(segments_array);
        domain.segments = segments.empty() ? nullptr : segments.data();
        domain.num_segments = segments.size();
    }

    if (!holes_object.is_none()) {
        auto holes_array = py::cast<py::array_t<double, py::array::c_style | py::array::forcecast>>(holes_object);
        if (holes_array.ndim() != 2 || holes_array.shape(1) != 2) {
            throw std::invalid_argument("holes must have shape (K, 2)");
        }
        holes = tm_points_from_array(holes_array);
        domain.holes = holes.empty() ? nullptr : holes.data();
        domain.num_holes = holes.size();
    }

    tm_options_init(&options);
    options.min_angle_deg = min_angle_deg;
    options.enable_refinement = enable_refinement ? 1 : 0;
    options.use_offcenters = use_offcenters ? 1 : 0;
    options.verbose = verbose ? 1 : 0;
    options.enable_acute_protection = acute_mode != 0 ? 1 : 0;
    options.acute_protection_mode = acute_mode == 1 ? TM_ACUTE_PROTECTION_SIMPLE : TM_ACUTE_PROTECTION_SHELL;
    options.max_refinement_steps = max_refinement_steps;
    options.max_protection_levels = max_protection_levels;
    if (!protect_angle_object.is_none()) {
        options.protect_angle_deg = py::cast<double>(protect_angle_object);
    }

    status = tm_generate(&domain, &options, &mesh, &error);
    if (status != TM_STATUS_OK) {
        throw std::runtime_error(error.message[0] != '\0' ? error.message : tm_status_string(status));
    }

    status = tm_analyze_mesh(&mesh, &summary, &error);
    if (status != TM_STATUS_OK) {
        tm_mesh_free(&mesh);
        throw std::runtime_error(error.message[0] != '\0' ? error.message : tm_status_string(status));
    }

    py::dict result;
    result["points"] = tm_points_to_numpy(mesh);
    result["triangles"] = tm_triangles_to_numpy(mesh);
    result["segments"] = tm_segments_to_numpy(mesh);
    result["summary"] = tm_summary_to_dict(summary);

    tm_mesh_free(&mesh);
    return result;
}

static py::dict tm_analyze_raw(
    const py::array_t<double, py::array::c_style | py::array::forcecast> &points_array,
    const py::array_t<std::uint32_t, py::array::c_style | py::array::forcecast> &triangles_array,
    py::object segments_object,
    std::size_t num_input_points,
    std::size_t num_segment_split_points,
    std::size_t num_triangle_split_points,
    std::size_t num_protected_corners,
    std::size_t num_exempt_triangles
)
{
    std::vector<tm_point> points;
    std::vector<tm_segment> segments;
    std::vector<std::uint32_t> triangles;
    tm_mesh mesh{};
    tm_quality_summary summary{};
    tm_error error{};
    tm_status status;

    if (points_array.ndim() != 2 || points_array.shape(1) != 2) {
        throw std::invalid_argument("points must have shape (N, 2)");
    }
    if (triangles_array.ndim() != 2 || triangles_array.shape(1) != 3) {
        throw std::invalid_argument("triangles must have shape (M, 3)");
    }

    points = tm_points_from_array(points_array);
    triangles.resize(static_cast<size_t>(triangles_array.shape(0)) * 3);
    {
        auto view = triangles_array.unchecked<2>();
        for (py::ssize_t i = 0; i < view.shape(0); ++i) {
            triangles[static_cast<size_t>(i) * 3 + 0] = view(i, 0);
            triangles[static_cast<size_t>(i) * 3 + 1] = view(i, 1);
            triangles[static_cast<size_t>(i) * 3 + 2] = view(i, 2);
        }
    }

    mesh.points = points.empty() ? nullptr : points.data();
    mesh.num_points = points.size();
    mesh.triangles = triangles.empty() ? nullptr : triangles.data();
    mesh.num_triangles = triangles.size() / 3;
    mesh.num_input_points = num_input_points;
    mesh.num_segment_split_points = num_segment_split_points;
    mesh.num_triangle_split_points = num_triangle_split_points;
    mesh.num_protected_corners = num_protected_corners;
    mesh.num_exempt_triangles = num_exempt_triangles;

    if (!segments_object.is_none()) {
        auto segments_array = py::cast<py::array_t<std::uint32_t, py::array::c_style | py::array::forcecast>>(segments_object);
        if (segments_array.ndim() != 2 || segments_array.shape(1) != 2) {
            throw std::invalid_argument("segments must have shape (M, 2)");
        }
        segments = tm_segments_from_array(segments_array);
        mesh.segments = segments.empty() ? nullptr : segments.data();
        mesh.num_segments = segments.size();
    }

    status = tm_analyze_mesh(&mesh, &summary, &error);
    if (status != TM_STATUS_OK) {
        throw std::runtime_error(error.message[0] != '\0' ? error.message : tm_status_string(status));
    }

    return tm_summary_to_dict(summary);
}

PYBIND11_MODULE(_core, m)
{
    m.doc() = "pybind11 bindings for the DTCC Mesher C library";
    m.attr("__version__") = DTCC_MESHER_VERSION_STRING;

    m.def(
        "_generate_raw",
        &tm_generate_raw,
        py::arg("points"),
        py::arg("segments") = py::none(),
        py::arg("holes") = py::none(),
        py::arg("min_angle_deg") = 20.0,
        py::arg("enable_refinement") = true,
        py::arg("use_offcenters") = false,
        py::arg("verbose") = false,
        py::arg("acute_mode") = 2,
        py::arg("protect_angle_deg") = py::none(),
        py::arg("max_refinement_steps") = 0,
        py::arg("max_protection_levels") = 6
    );
    m.def(
        "_analyze_raw",
        &tm_analyze_raw,
        py::arg("points"),
        py::arg("triangles"),
        py::arg("segments") = py::none(),
        py::arg("num_input_points") = 0,
        py::arg("num_segment_split_points") = 0,
        py::arg("num_triangle_split_points") = 0,
        py::arg("num_protected_corners") = 0,
        py::arg("num_exempt_triangles") = 0
    );
}
