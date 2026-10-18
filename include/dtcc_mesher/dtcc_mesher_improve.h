#ifndef DTCC_MESHER_IMPROVE_H
#define DTCC_MESHER_IMPROVE_H
#include "dtcc_mesher.h"
#ifdef __cplusplus
extern "C" {
#endif

/* Shape optimization is stateless and may produce non-Delaunay triangles.
 * Input vertices/indices are immutable. max_added_points is exactly 0 or 1.
 * Only listed input segment indices may split, including internal constraints.
 * Every boundary edge must be a segment. Input must be a conforming oriented
 * planar triangle mesh; unsupported/nonmanifold/overlapping inputs are errors.
 * Initial scope: at most 4096 points and 4096 triangles per planar domain.
 */
typedef enum {
    DTCC_MESHER_INSERT_NONE = 0,
    DTCC_MESHER_INSERT_INTERIOR = 1,
    DTCC_MESHER_INSERT_SEGMENT = 2
} dtcc_mesher_insertion_kind;

typedef struct {
    double penalty;                  /* sum((1-q)^2), q=4*sqrt(3)*area/sum(edge^2) */
    double minimum_q;
    size_t poor_5_count;
    double poor_5_area;
    size_t poor_10_count;
    double poor_10_area;
} dtcc_mesher_shape_score;

#define DTCC_MESHER_SEGMENT_TRIAL_LIMIT 4
typedef struct {
    size_t parent_segment;
    double segment_parameter;
} dtcc_mesher_segment_trial;

typedef struct {
    dtcc_mesher_mesh mesh;
    size_t *segment_parents;         /* one input segment index per output segment */
    dtcc_mesher_insertion_kind insertion_kind;
    size_t inserted_point;           /* SIZE_MAX when no point was inserted */
    size_t parent_segment;           /* SIZE_MAX except for segment insertion */
    double segment_parameter;        /* input segment a + t*(b-a), 0 < t < 1 */
    size_t flips;
    size_t candidate_evaluations;
    dtcc_mesher_shape_score before;
    dtcc_mesher_shape_score after;
    dtcc_mesher_quality_summary quality_before;
    dtcc_mesher_quality_summary quality_after;
    char termination[64];
    /* Valid beneficial one-point segment trials, ranked by gain/triangle.
     * Includes the winner when it is a segment trial. No additional meshes
     * or allocations: callers may reconcile one parameter with neighbours.
     */
    size_t num_segment_trials;
    dtcc_mesher_segment_trial segment_trials[DTCC_MESHER_SEGMENT_TRIAL_LIMIT];
} dtcc_mesher_improvement;

/* Zero-initialize an empty result. Never alias input and output storage.
 * The operation borrows input. On success the result owns all returned arrays;
 * on failure it is empty. Free even unchanged results with improvement_free.
 */
DTCC_MESHER_EXPORT void dtcc_mesher_improvement_init(dtcc_mesher_improvement *result);
DTCC_MESHER_EXPORT void dtcc_mesher_improvement_free(dtcc_mesher_improvement *result);
DTCC_MESHER_EXPORT dtcc_mesher_status dtcc_mesher_improve(
    const dtcc_mesher_mesh *input,
    unsigned int max_added_points,
    int allow_interior,
    const size_t *splittable_segments,
    size_t num_splittable_segments,
    dtcc_mesher_improvement *result,
    dtcc_mesher_error *error);
#ifdef __cplusplus
}
#endif
#endif
