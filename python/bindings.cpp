#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "dtcc_mesher/dtcc_mesher.h"
#include "dtcc_mesher/dtcc_mesher_version.h"
extern "C" {
#include "cdt.h"
#include "io_pslg.h"
}

namespace py = pybind11;

static py::dict dtcc_mesher_summary_to_dict(const dtcc_mesher_quality_summary &summary);

static std::vector<dtcc_mesher_point> dtcc_mesher_points_from_array(const py::array_t<double, py::array::c_style | py::array::forcecast> &array)
{
    auto view = array.unchecked<2>();
    std::vector<dtcc_mesher_point> points(static_cast<size_t>(view.shape(0)));

    for (py::ssize_t i = 0; i < view.shape(0); ++i) {
        points[static_cast<size_t>(i)].x = view(i, 0);
        points[static_cast<size_t>(i)].y = view(i, 1);
    }

    return points;
}

static std::vector<dtcc_mesher_segment> dtcc_mesher_segments_from_array(const py::array_t<std::uint32_t, py::array::c_style | py::array::forcecast> &array)
{
    auto view = array.unchecked<2>();
    std::vector<dtcc_mesher_segment> segments(static_cast<size_t>(view.shape(0)));

    for (py::ssize_t i = 0; i < view.shape(0); ++i) {
        segments[static_cast<size_t>(i)].a = view(i, 0);
        segments[static_cast<size_t>(i)].b = view(i, 1);
    }

    return segments;
}

static std::vector<TMPoint> dtcc_mesher_internal_points_from_array(const py::array_t<double, py::array::c_style | py::array::forcecast> &array)
{
    auto view = array.unchecked<2>();
    std::vector<TMPoint> points(static_cast<size_t>(view.shape(0)));

    for (py::ssize_t i = 0; i < view.shape(0); ++i) {
        points[static_cast<size_t>(i)].xy[0] = view(i, 0);
        points[static_cast<size_t>(i)].xy[1] = view(i, 1);
        points[static_cast<size_t>(i)].original_index = static_cast<int>(i);
        points[static_cast<size_t>(i)].kind = TM_VERTEX_INPUT;
        points[static_cast<size_t>(i)].incident_triangle = -1;
        points[static_cast<size_t>(i)].protection_apex = -1;
        points[static_cast<size_t>(i)].protection_side = TM_PROTECTION_SIDE_NONE;
        points[static_cast<size_t>(i)].protection_level = 0;
    }

    return points;
}

static std::vector<TMSegment> dtcc_mesher_internal_segments_from_array(const py::array_t<std::uint32_t, py::array::c_style | py::array::forcecast> &array)
{
    auto view = array.unchecked<2>();
    std::vector<TMSegment> segments(static_cast<size_t>(view.shape(0)));

    for (py::ssize_t i = 0; i < view.shape(0); ++i) {
        segments[static_cast<size_t>(i)].v[0] = static_cast<int>(view(i, 0));
        segments[static_cast<size_t>(i)].v[1] = static_cast<int>(view(i, 1));
        segments[static_cast<size_t>(i)].original_index = static_cast<int>(i);
        segments[static_cast<size_t>(i)].live = 1;
        segments[static_cast<size_t>(i)].is_protected = 0;
        segments[static_cast<size_t>(i)].protected_apex = -1;
    }

    return segments;
}

static std::vector<TMRegion> dtcc_mesher_internal_regions_from_arrays(
    const py::array_t<double, py::array::c_style | py::array::forcecast> &points_array,
    const py::array_t<std::int32_t, py::array::c_style | py::array::forcecast> &markers_array
)
{
    auto points_view = points_array.unchecked<2>();
    auto markers_view = markers_array.unchecked<1>();
    std::vector<TMRegion> regions(static_cast<size_t>(points_view.shape(0)));

    for (py::ssize_t i = 0; i < points_view.shape(0); ++i) {
        regions[static_cast<size_t>(i)].xy[0] = points_view(i, 0);
        regions[static_cast<size_t>(i)].xy[1] = points_view(i, 1);
        regions[static_cast<size_t>(i)].marker = markers_view(i);
    }

    return regions;
}

static py::array_t<double> dtcc_mesher_points_to_numpy(const dtcc_mesher_mesh &mesh)
{
    py::array_t<double> array({static_cast<py::ssize_t>(mesh.num_points), static_cast<py::ssize_t>(2)});
    auto view = array.mutable_unchecked<2>();

    for (size_t i = 0; i < mesh.num_points; ++i) {
        view(static_cast<py::ssize_t>(i), 0) = mesh.points[i].x;
        view(static_cast<py::ssize_t>(i), 1) = mesh.points[i].y;
    }

    return array;
}

static py::array_t<std::uint32_t> dtcc_mesher_triangles_to_numpy(const dtcc_mesher_mesh &mesh)
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

static py::array_t<std::uint32_t> dtcc_mesher_segments_to_numpy(const dtcc_mesher_mesh &mesh)
{
    py::array_t<std::uint32_t> array({static_cast<py::ssize_t>(mesh.num_segments), static_cast<py::ssize_t>(2)});
    auto view = array.mutable_unchecked<2>();

    for (size_t i = 0; i < mesh.num_segments; ++i) {
        view(static_cast<py::ssize_t>(i), 0) = mesh.segments[i].a;
        view(static_cast<py::ssize_t>(i), 1) = mesh.segments[i].b;
    }

    return array;
}

static py::array_t<std::int32_t> dtcc_mesher_markers_to_numpy(const int *markers, size_t count)
{
    py::array_t<std::int32_t> array(static_cast<py::ssize_t>(count));
    auto view = array.mutable_unchecked<1>();

    for (size_t i = 0; i < count; ++i) {
        view(static_cast<py::ssize_t>(i)) = static_cast<std::int32_t>(markers[i]);
    }

    return array;
}

static void dtcc_mesher_options_from_python(
    dtcc_mesher_options *options,
    double min_angle_deg,
    py::object max_area_object,
    py::object max_edge_length_object,
    bool enable_refinement,
    bool verbose,
    int acute_mode,
    py::object protect_angle_object,
    std::size_t max_refinement_steps,
    std::size_t max_protection_levels
)
{
    dtcc_mesher_options_init(options);
    options->min_angle_deg = min_angle_deg;
    if (!max_area_object.is_none()) {
        options->max_area = py::cast<double>(max_area_object);
    }
    if (!max_edge_length_object.is_none()) {
        options->max_edge_length = py::cast<double>(max_edge_length_object);
    }
    options->enable_refinement = enable_refinement ? 1 : 0;
    options->verbose = verbose ? 1 : 0;
    options->enable_acute_protection = acute_mode != 0 ? 1 : 0;
    options->acute_protection_mode = acute_mode == 1 ? DTCC_MESHER_ACUTE_PROTECTION_SIMPLE : DTCC_MESHER_ACUTE_PROTECTION_SHELL;
    options->max_refinement_steps = max_refinement_steps;
    options->max_protection_levels = max_protection_levels;
    if (!protect_angle_object.is_none()) {
        options->protect_angle_deg = py::cast<double>(protect_angle_object);
    }
}

static TMBuildOptions dtcc_mesher_build_options_from_public(const dtcc_mesher_options &options)
{
    TMBuildOptions build_options{};

    build_options.verbose = options.verbose;
    build_options.refine = options.enable_refinement || options.max_area > 0.0;
    build_options.protect_acute_corners = options.enable_acute_protection;
    build_options.acute_mode = options.acute_protection_mode == DTCC_MESHER_ACUTE_PROTECTION_SIMPLE ?
        TM_ACUTE_MODE_SIMPLE :
        TM_ACUTE_MODE_SHELL;
    build_options.min_angle_deg = options.min_angle_deg;
    build_options.max_area = options.max_area;
    build_options.max_edge_length = options.max_edge_length;
    build_options.protect_angle_deg = options.protect_angle_deg;
    build_options.max_refinement_steps = options.max_refinement_steps;
    build_options.max_protection_levels = options.max_protection_levels;
    return build_options;
}

static void dtcc_mesher_internal_mesh_metadata(const TMMesh &mesh, dtcc_mesher_mesh *public_mesh)
{
    size_t point_index;
    size_t triangle_index;

    public_mesh->num_input_points = 0;
    public_mesh->num_segment_split_points = 0;
    public_mesh->num_triangle_split_points = 0;
    public_mesh->num_protected_corners = tm_count_protected_corners(&mesh);
    public_mesh->num_exempt_triangles = 0;

    for (point_index = 0; point_index < mesh.point_count; ++point_index) {
        switch (mesh.points[point_index].kind) {
            case TM_VERTEX_INPUT:
                public_mesh->num_input_points += 1;
                break;
            case TM_VERTEX_SEGMENT_SPLIT:
                public_mesh->num_segment_split_points += 1;
                break;
            case TM_VERTEX_TRIANGLE_SPLIT:
                public_mesh->num_triangle_split_points += 1;
                break;
            case TM_VERTEX_SUPER:
                break;
        }
    }

    for (triangle_index = 0; triangle_index < mesh.triangle_count; ++triangle_index) {
        if (tm_triangle_is_quality_exempt(&mesh, static_cast<int>(triangle_index))) {
            public_mesh->num_exempt_triangles += 1;
        }
    }
}

static py::dict dtcc_mesher_internal_mesh_result(const TMMesh &mesh, const int *triangle_markers)
{
    std::vector<dtcc_mesher_point> points(mesh.point_count);
    std::vector<dtcc_mesher_segment> segments;
    std::vector<std::uint32_t> triangles(mesh.triangle_count * 3);
    dtcc_mesher_mesh public_mesh{};
    dtcc_mesher_quality_summary summary{};
    dtcc_mesher_error error{};
    dtcc_mesher_status status;
    py::array_t<double> points_array({static_cast<py::ssize_t>(points.size()), static_cast<py::ssize_t>(2)});
    py::array_t<std::uint32_t> triangles_array({static_cast<py::ssize_t>(mesh.triangle_count), static_cast<py::ssize_t>(3)});
    py::array_t<std::uint32_t> segments_array({static_cast<py::ssize_t>(segments.size()), static_cast<py::ssize_t>(2)});
    size_t i;

    for (i = 0; i < mesh.point_count; ++i) {
        points[i].x = mesh.points[i].xy[0];
        points[i].y = mesh.points[i].xy[1];
    }

    for (i = 0; i < mesh.triangle_count; ++i) {
        triangles[i * 3 + 0] = static_cast<std::uint32_t>(mesh.triangles[i].v[0]);
        triangles[i * 3 + 1] = static_cast<std::uint32_t>(mesh.triangles[i].v[1]);
        triangles[i * 3 + 2] = static_cast<std::uint32_t>(mesh.triangles[i].v[2]);
    }

    for (i = 0; i < mesh.segment_count; ++i) {
        if (!mesh.segments[i].live) {
            continue;
        }
        segments.push_back(
            dtcc_mesher_segment{
                static_cast<std::uint32_t>(mesh.segments[i].v[0]),
                static_cast<std::uint32_t>(mesh.segments[i].v[1]),
            }
        );
    }

    public_mesh.points = points.empty() ? nullptr : points.data();
    public_mesh.num_points = points.size();
    public_mesh.triangles = triangles.empty() ? nullptr : triangles.data();
    public_mesh.num_triangles = triangles.size() / 3;
    public_mesh.segments = segments.empty() ? nullptr : segments.data();
    public_mesh.num_segments = segments.size();
    dtcc_mesher_internal_mesh_metadata(mesh, &public_mesh);

    status = dtcc_mesher_analyze_mesh(&public_mesh, &summary, &error);
    if (status != DTCC_MESHER_STATUS_OK) {
        throw std::runtime_error(error.message[0] != '\0' ? error.message : dtcc_mesher_status_string(status));
    }

    py::dict result;
    {
        auto view = points_array.mutable_unchecked<2>();
        for (i = 0; i < points.size(); ++i) {
            view(static_cast<py::ssize_t>(i), 0) = points[i].x;
            view(static_cast<py::ssize_t>(i), 1) = points[i].y;
        }
    }
    {
        auto view = triangles_array.mutable_unchecked<2>();
        for (i = 0; i < mesh.triangle_count; ++i) {
            view(static_cast<py::ssize_t>(i), 0) = triangles[i * 3 + 0];
            view(static_cast<py::ssize_t>(i), 1) = triangles[i * 3 + 1];
            view(static_cast<py::ssize_t>(i), 2) = triangles[i * 3 + 2];
        }
    }
    {
        auto view = segments_array.mutable_unchecked<2>();
        for (i = 0; i < segments.size(); ++i) {
            view(static_cast<py::ssize_t>(i), 0) = segments[i].a;
            view(static_cast<py::ssize_t>(i), 1) = segments[i].b;
        }
    }
    result["points"] = points_array;
    result["triangles"] = triangles_array;
    result["segments"] = segments_array;
    if (triangle_markers != nullptr) {
        result["markers"] = dtcc_mesher_markers_to_numpy(triangle_markers, mesh.triangle_count);
    }
    result["summary"] = dtcc_mesher_summary_to_dict(summary);
    return result;
}

static py::dict dtcc_mesher_summary_to_dict(const dtcc_mesher_quality_summary &summary)
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

static py::dict dtcc_mesher_generate_raw(
    const py::array_t<double, py::array::c_style | py::array::forcecast> &points_array,
    py::object segments_object,
    py::object holes_object,
    double min_angle_deg,
    py::object max_area_object,
    py::object max_edge_length_object,
    bool enable_refinement,
    bool verbose,
    int acute_mode,
    py::object protect_angle_object,
    std::size_t max_refinement_steps,
    std::size_t max_protection_levels
)
{
    std::vector<dtcc_mesher_point> points;
    std::vector<dtcc_mesher_segment> segments;
    std::vector<dtcc_mesher_point> holes;
    dtcc_mesher_domain domain{};
    dtcc_mesher_options options;
    dtcc_mesher_mesh mesh{};
    dtcc_mesher_quality_summary summary{};
    dtcc_mesher_error error{};
    dtcc_mesher_status status;

    if (points_array.ndim() != 2 || points_array.shape(1) != 2) {
        throw std::invalid_argument("points must have shape (N, 2)");
    }

    points = dtcc_mesher_points_from_array(points_array);
    domain.points = points.empty() ? nullptr : points.data();
    domain.num_points = points.size();

    if (!segments_object.is_none()) {
        auto segments_array = py::cast<py::array_t<std::uint32_t, py::array::c_style | py::array::forcecast>>(segments_object);
        if (segments_array.ndim() != 2 || segments_array.shape(1) != 2) {
            throw std::invalid_argument("segments must have shape (M, 2)");
        }
        segments = dtcc_mesher_segments_from_array(segments_array);
        domain.segments = segments.empty() ? nullptr : segments.data();
        domain.num_segments = segments.size();
    }

    if (!holes_object.is_none()) {
        auto holes_array = py::cast<py::array_t<double, py::array::c_style | py::array::forcecast>>(holes_object);
        if (holes_array.ndim() != 2 || holes_array.shape(1) != 2) {
            throw std::invalid_argument("holes must have shape (K, 2)");
        }
        holes = dtcc_mesher_points_from_array(holes_array);
        domain.holes = holes.empty() ? nullptr : holes.data();
        domain.num_holes = holes.size();
    }

    dtcc_mesher_options_from_python(
        &options,
        min_angle_deg,
        max_area_object,
        max_edge_length_object,
        enable_refinement,
        verbose,
        acute_mode,
        protect_angle_object,
        max_refinement_steps,
        max_protection_levels
    );

    status = dtcc_mesher_generate(&domain, &options, &mesh, &error);
    if (status != DTCC_MESHER_STATUS_OK) {
        throw std::runtime_error(error.message[0] != '\0' ? error.message : dtcc_mesher_status_string(status));
    }

    status = dtcc_mesher_analyze_mesh(&mesh, &summary, &error);
    if (status != DTCC_MESHER_STATUS_OK) {
        dtcc_mesher_mesh_free(&mesh);
        throw std::runtime_error(error.message[0] != '\0' ? error.message : dtcc_mesher_status_string(status));
    }

    py::dict result;
    result["points"] = dtcc_mesher_points_to_numpy(mesh);
    result["triangles"] = dtcc_mesher_triangles_to_numpy(mesh);
    result["segments"] = dtcc_mesher_segments_to_numpy(mesh);
    result["summary"] = dtcc_mesher_summary_to_dict(summary);

    dtcc_mesher_mesh_free(&mesh);
    return result;
}

static py::dict dtcc_mesher_generate_coverage_raw(
    const py::array_t<double, py::array::c_style | py::array::forcecast> &points_array,
    const py::array_t<std::uint32_t, py::array::c_style | py::array::forcecast> &segments_array,
    const py::array_t<double, py::array::c_style | py::array::forcecast> &region_points_array,
    const py::array_t<std::int32_t, py::array::c_style | py::array::forcecast> &region_markers_array,
    double min_angle_deg,
    py::object max_area_object,
    py::object max_edge_length_object,
    bool enable_refinement,
    bool verbose,
    int acute_mode,
    py::object protect_angle_object,
    std::size_t max_refinement_steps,
    std::size_t max_protection_levels
)
{
    std::vector<TMPoint> points;
    std::vector<TMSegment> segments;
    std::vector<TMRegion> regions;
    TMPSLG pslg{};
    TMMesh mesh{};
    dtcc_mesher_options public_options{};
    TMBuildOptions build_options{};
    int *triangle_markers = nullptr;
    TMStatus status;

    if (points_array.ndim() != 2 || points_array.shape(1) != 2) {
        throw std::invalid_argument("points must have shape (N, 2)");
    }
    if (segments_array.ndim() != 2 || segments_array.shape(1) != 2) {
        throw std::invalid_argument("segments must have shape (M, 2)");
    }
    if (region_points_array.ndim() != 2 || region_points_array.shape(1) != 2) {
        throw std::invalid_argument("region_points must have shape (R, 2)");
    }
    if (region_markers_array.ndim() != 1) {
        throw std::invalid_argument("region_markers must have shape (R,)");
    }
    if (region_points_array.shape(0) != region_markers_array.shape(0)) {
        throw std::invalid_argument("region_points and region_markers must have the same length");
    }

    points = dtcc_mesher_internal_points_from_array(points_array);
    segments = dtcc_mesher_internal_segments_from_array(segments_array);
    regions = dtcc_mesher_internal_regions_from_arrays(region_points_array, region_markers_array);
    pslg.points = points.empty() ? nullptr : points.data();
    pslg.point_count = points.size();
    pslg.segments = segments.empty() ? nullptr : segments.data();
    pslg.segment_count = segments.size();
    pslg.holes = nullptr;
    pslg.hole_count = 0;

    dtcc_mesher_options_from_python(
        &public_options,
        min_angle_deg,
        max_area_object,
        max_edge_length_object,
        enable_refinement,
        verbose,
        acute_mode,
        protect_angle_object,
        max_refinement_steps,
        max_protection_levels
    );
    build_options = dtcc_mesher_build_options_from_public(public_options);

    if (public_options.min_angle_deg <= 0.0 || public_options.min_angle_deg >= 90.0) {
        throw std::invalid_argument("min_angle_deg must be in the range (0, 90)");
    }
    if (public_options.max_area < 0.0) {
        throw std::invalid_argument("max_area must be non-negative");
    }
    if (public_options.max_edge_length < 0.0) {
        throw std::invalid_argument("max_edge_length must be non-negative");
    }
    if (public_options.protect_angle_deg < 0.0 || public_options.protect_angle_deg >= 180.0) {
        throw std::invalid_argument("protect_angle_deg must be in the range [0, 180)");
    }

    status = tm_build_coverage_mesh(
        &pslg,
        regions.empty() ? nullptr : regions.data(),
        regions.size(),
        &build_options,
        &mesh,
        &triangle_markers
    );
    if (status != TM_OK) {
        const char *detail = tm_last_pslg_error_detail();
        tm_free_mesh(&mesh);
        free(triangle_markers);
        throw std::runtime_error(detail != nullptr && detail[0] != '\0' ? detail : tm_internal_status_string(status));
    }

    try {
        py::dict result = dtcc_mesher_internal_mesh_result(mesh, triangle_markers);
        tm_free_mesh(&mesh);
        free(triangle_markers);
        return result;
    } catch (...) {
        tm_free_mesh(&mesh);
        free(triangle_markers);
        throw;
    }
}

static void dtcc_mesher_validate_segment_graph_raw(
    const py::array_t<double, py::array::c_style | py::array::forcecast> &points_array,
    const py::array_t<std::uint32_t, py::array::c_style | py::array::forcecast> &segments_array
)
{
    std::vector<TMPoint> points;
    std::vector<TMSegment> segments;
    TMStatus status;

    if (points_array.ndim() != 2 || points_array.shape(1) != 2) {
        throw std::invalid_argument("points must have shape (N, 2)");
    }
    if (segments_array.ndim() != 2 || segments_array.shape(1) != 2) {
        throw std::invalid_argument("segments must have shape (M, 2)");
    }

    points = dtcc_mesher_internal_points_from_array(points_array);
    segments = dtcc_mesher_internal_segments_from_array(segments_array);
    status = tm_validate_segment_graph(
        points.empty() ? nullptr : points.data(),
        points.size(),
        segments.empty() ? nullptr : segments.data(),
        segments.size()
    );
    if (status != TM_OK) {
        const char *detail = tm_last_pslg_error_detail();
        throw std::runtime_error(
            detail != nullptr && detail[0] != '\0'
                ? detail
                : tm_internal_status_string(status)
        );
    }
}

static py::dict dtcc_mesher_analyze_raw(
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
    std::vector<dtcc_mesher_point> points;
    std::vector<dtcc_mesher_segment> segments;
    std::vector<std::uint32_t> triangles;
    dtcc_mesher_mesh mesh{};
    dtcc_mesher_quality_summary summary{};
    dtcc_mesher_error error{};
    dtcc_mesher_status status;

    if (points_array.ndim() != 2 || points_array.shape(1) != 2) {
        throw std::invalid_argument("points must have shape (N, 2)");
    }
    if (triangles_array.ndim() != 2 || triangles_array.shape(1) != 3) {
        throw std::invalid_argument("triangles must have shape (M, 3)");
    }

    points = dtcc_mesher_points_from_array(points_array);
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
        segments = dtcc_mesher_segments_from_array(segments_array);
        mesh.segments = segments.empty() ? nullptr : segments.data();
        mesh.num_segments = segments.size();
    }

    status = dtcc_mesher_analyze_mesh(&mesh, &summary, &error);
    if (status != DTCC_MESHER_STATUS_OK) {
        throw std::runtime_error(error.message[0] != '\0' ? error.message : dtcc_mesher_status_string(status));
    }

    return dtcc_mesher_summary_to_dict(summary);
}

PYBIND11_MODULE(_core, m)
{
    m.doc() = "pybind11 bindings for the DTCC Mesher C library";
    m.attr("__version__") = DTCC_MESHER_VERSION_STRING;

    m.def(
        "_generate_raw",
        &dtcc_mesher_generate_raw,
        py::arg("points"),
        py::arg("segments") = py::none(),
        py::arg("holes") = py::none(),
        py::arg("min_angle_deg") = 20.0,
        py::arg("max_area") = py::none(),
        py::arg("max_edge_length") = py::none(),
        py::arg("enable_refinement") = true,
        py::arg("verbose") = false,
        py::arg("acute_mode") = 2,
        py::arg("protect_angle_deg") = py::none(),
        py::arg("max_refinement_steps") = 0,
        py::arg("max_protection_levels") = 6
    );
    m.def(
        "_generate_coverage_raw",
        &dtcc_mesher_generate_coverage_raw,
        py::arg("points"),
        py::arg("segments"),
        py::arg("region_points"),
        py::arg("region_markers"),
        py::arg("min_angle_deg") = 20.0,
        py::arg("max_area") = py::none(),
        py::arg("max_edge_length") = py::none(),
        py::arg("enable_refinement") = true,
        py::arg("verbose") = false,
        py::arg("acute_mode") = 2,
        py::arg("protect_angle_deg") = py::none(),
        py::arg("max_refinement_steps") = 0,
        py::arg("max_protection_levels") = 6
    );
    m.def(
        "_validate_segment_graph_raw",
        &dtcc_mesher_validate_segment_graph_raw,
        py::arg("points"),
        py::arg("segments")
    );
    m.def(
        "_analyze_raw",
        &dtcc_mesher_analyze_raw,
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
