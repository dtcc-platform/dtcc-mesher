#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dtcc_mesher/dtcc_mesher.h"
#include "dtcc_mesher/dtcc_mesher_io.h"
#include "dtcc_mesher/dtcc_mesher_version.h"

static int tm_has_suffix(const char *path, const char *suffix)
{
    size_t path_len = strlen(path);
    size_t suffix_len = strlen(suffix);

    return path_len >= suffix_len && strcmp(path + path_len - suffix_len, suffix) == 0;
}

static char *tm_with_suffix(const char *base, const char *suffix)
{
    char *path;
    size_t base_len = strlen(base);
    size_t suffix_len = strlen(suffix);

    path = (char *) malloc(base_len + suffix_len + 1);
    if (path == NULL) {
        return NULL;
    }

    memcpy(path, base, base_len);
    memcpy(path + base_len, suffix, suffix_len + 1);
    return path;
}

static int tm_parse_double_arg(const char *text, double *out_value)
{
    char *end = NULL;
    double value = strtod(text, &end);

    if (end == text || *end != '\0') {
        return 0;
    }

    *out_value = value;
    return 1;
}

static int tm_parse_size_arg(const char *text, size_t *out_value)
{
    char *end = NULL;
    unsigned long long value = strtoull(text, &end, 10);

    if (end == text || *end != '\0') {
        return 0;
    }

    *out_value = (size_t) value;
    return 1;
}

static void tm_print_usage(FILE *stream)
{
    fprintf(
        stream,
        "usage: dtcc_mesher [-h|--help] [--version] [-v|--verbose] [--refine|--no-refine] [--off-centers|--no-off-centers] [--acute-protection|--no-acute-protection] [--simple-acute-protection|--shell-acute-protection] [--min-angle deg] [--protect-angle deg] [--max-refine-steps n] [--max-protection-levels n] input.(pts|pslg) outbase\n"
    );
}

static int tm_print_error(const char *context, const char *path, tm_status status, const tm_error *error)
{
    const char *detail = tm_status_string(status);

    if (error != NULL && error->message[0] != '\0') {
        detail = error->message;
    }

    if (path != NULL) {
        fprintf(stderr, "%s %s: %s\n", context, path, detail);
    } else {
        fprintf(stderr, "%s: %s\n", context, detail);
    }

    return EXIT_FAILURE;
}

int main(int argc, char **argv)
{
    const char *input_path = NULL;
    const char *out_base = NULL;
    tm_options options;
    tm_domain domain;
    tm_mesh mesh;
    tm_error error;
    char *tri_path = NULL;
    char *svg_path = NULL;
    char *metrics_path = NULL;
    char *summary_path = NULL;
    tm_status status;
    int exit_code = EXIT_FAILURE;
    int argi = 1;

    tm_options_init(&options);
    memset(&domain, 0, sizeof(domain));
    memset(&mesh, 0, sizeof(mesh));
    memset(&error, 0, sizeof(error));

    while (argi < argc && argv[argi][0] == '-') {
        if (strcmp(argv[argi], "-h") == 0 || strcmp(argv[argi], "--help") == 0) {
            tm_print_usage(stdout);
            return EXIT_SUCCESS;
        }
        if (strcmp(argv[argi], "--version") == 0) {
            printf("dtcc_mesher %s\n", DTCC_MESHER_VERSION_STRING);
            return EXIT_SUCCESS;
        }
        if (strcmp(argv[argi], "-v") == 0 || strcmp(argv[argi], "--verbose") == 0) {
            options.verbose = 1;
        } else if (strcmp(argv[argi], "--refine") == 0) {
            options.enable_refinement = 1;
        } else if (strcmp(argv[argi], "--no-refine") == 0) {
            options.enable_refinement = 0;
        } else if (strcmp(argv[argi], "--off-centers") == 0) {
            options.use_offcenters = 1;
        } else if (strcmp(argv[argi], "--no-off-centers") == 0) {
            options.use_offcenters = 0;
        } else if (strcmp(argv[argi], "--acute-protection") == 0) {
            options.enable_acute_protection = 1;
        } else if (strcmp(argv[argi], "--no-acute-protection") == 0) {
            options.enable_acute_protection = 0;
            options.acute_protection_mode = TM_ACUTE_PROTECTION_NONE;
        } else if (strcmp(argv[argi], "--simple-acute-protection") == 0) {
            options.enable_acute_protection = 1;
            options.acute_protection_mode = TM_ACUTE_PROTECTION_SIMPLE;
        } else if (strcmp(argv[argi], "--shell-acute-protection") == 0) {
            options.enable_acute_protection = 1;
            options.acute_protection_mode = TM_ACUTE_PROTECTION_SHELL;
        } else if (strcmp(argv[argi], "--min-angle") == 0) {
            argi += 1;
            if (argi >= argc || !tm_parse_double_arg(argv[argi], &options.min_angle_deg)) {
                tm_print_usage(stderr);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[argi], "--protect-angle") == 0) {
            argi += 1;
            if (argi >= argc || !tm_parse_double_arg(argv[argi], &options.protect_angle_deg)) {
                tm_print_usage(stderr);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[argi], "--max-refine-steps") == 0) {
            argi += 1;
            if (argi >= argc || !tm_parse_size_arg(argv[argi], &options.max_refinement_steps)) {
                tm_print_usage(stderr);
                return EXIT_FAILURE;
            }
        } else if (strcmp(argv[argi], "--max-protection-levels") == 0) {
            argi += 1;
            if (argi >= argc || !tm_parse_size_arg(argv[argi], &options.max_protection_levels)) {
                tm_print_usage(stderr);
                return EXIT_FAILURE;
            }
        } else {
            tm_print_usage(stderr);
            return EXIT_FAILURE;
        }
        argi += 1;
    }

    if (argc - argi != 2) {
        tm_print_usage(stderr);
        return EXIT_FAILURE;
    }

    input_path = argv[argi];
    out_base = argv[argi + 1];

    if (!tm_has_suffix(input_path, ".pts") && !tm_has_suffix(input_path, ".pslg")) {
        fprintf(stderr, "unsupported input type for %s\n", input_path);
        return EXIT_FAILURE;
    }

    status = tm_read_domain_file(input_path, &domain, &error);
    if (status != TM_STATUS_OK) {
        return tm_print_error("failed to read", input_path, status, &error);
    }

    status = tm_generate(&domain, &options, &mesh, &error);
    if (status != TM_STATUS_OK) {
        tm_domain_free(&domain);
        return tm_print_error("failed to mesh", input_path, status, &error);
    }

    tri_path = tm_with_suffix(out_base, ".tri");
    svg_path = tm_with_suffix(out_base, ".svg");
    metrics_path = tm_with_suffix(out_base, ".metrics.csv");
    summary_path = tm_with_suffix(out_base, ".summary.txt");
    if (tri_path == NULL || svg_path == NULL || metrics_path == NULL || summary_path == NULL) {
        fprintf(stderr, "failed to allocate output paths\n");
        goto cleanup;
    }

    status = tm_write_triangles(&mesh, tri_path, &error);
    if (status != TM_STATUS_OK) {
        exit_code = tm_print_error("failed to write", tri_path, status, &error);
        goto cleanup;
    }

    status = tm_write_svg(&mesh, svg_path, &error);
    if (status != TM_STATUS_OK) {
        exit_code = tm_print_error("failed to write", svg_path, status, &error);
        goto cleanup;
    }

    status = tm_write_quality_csv(&mesh, metrics_path, &error);
    if (status != TM_STATUS_OK) {
        exit_code = tm_print_error("failed to write", metrics_path, status, &error);
        goto cleanup;
    }

    status = tm_write_quality_summary(&mesh, summary_path, &error);
    if (status != TM_STATUS_OK) {
        exit_code = tm_print_error("failed to write", summary_path, status, &error);
        goto cleanup;
    }

    exit_code = EXIT_SUCCESS;

cleanup:
    free(tri_path);
    free(svg_path);
    free(metrics_path);
    free(summary_path);
    tm_mesh_free(&mesh);
    tm_domain_free(&domain);
    return exit_code;
}
