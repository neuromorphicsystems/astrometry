#include "astrometry_extension_utilities.h"

typedef struct astrometry_extension_solver_t {
    PyObject_HEAD pl* indexes;
} astrometry_extension_solver_t;

static void astrometry_extension_solver_dealloc(PyObject* self) {
    astrometry_extension_solver_t* current = (astrometry_extension_solver_t*)self;
    if (current->indexes) {
        pl_remove_all(current->indexes);
        current->indexes = NULL;
    }
    Py_TYPE(self)->tp_free(self);
}

static PyObject* astrometry_extension_solver_new(PyTypeObject* type, PyObject* args, PyObject* kwds) {
    return type->tp_alloc(type, 0);
}

static int astrometry_extension_solver_init(PyObject* self, PyObject* args, PyObject* kwds) {
    PyObject* paths;
    if (!PyArg_ParseTuple(args, "O", &paths)) {
        return -1;
    }
    if (!PyList_Check(paths)) {
        PyErr_SetString(PyExc_TypeError, "paths must be a list");
        return -1;
    }
    if (PyList_GET_SIZE(paths) == 0) {
        PyErr_SetString(PyExc_TypeError, "paths cannot be empty");
        return -1;
    }
    astrometry_extension_solver_t* current = (astrometry_extension_solver_t*)self;
    current->indexes = pl_new(PyList_GET_SIZE(paths));
    for (Py_ssize_t path_index = 0; path_index < PyList_GET_SIZE(paths); ++path_index) {
        PyObject* path = PyList_GET_ITEM(paths, path_index);
        if (!PyUnicode_Check(path)) {
            PyErr_SetString(PyExc_TypeError, "all the items in paths must be strings");
            break;
        }
        const char* filename = (const char*)PyUnicode_DATA(path);
        anqfits_t* fits = anqfits_open(filename);
        if (fits == NULL) {
            PyErr_Format(PyExc_TypeError, "loading \"%s\" failed", filename);
            break;
        }
        index_t* index = calloc(1, sizeof(index_t));
        index->fits = fits;
        index->indexfn = (char*)filename;
        if (index_reload(index) != 0) {
            anqfits_close(index->fits);
            free(index);
            PyErr_Format(PyExc_TypeError, "loading \"%s\" failed", filename);
            break;
        }
        index->indexfn = strdup(index->indexfn);
        index->indexname = strdup(quadfile_get_filename(index->quads));
        index->index_scale_upper = quadfile_get_index_scale_upper_arcsec(index->quads);
        index->index_scale_lower = quadfile_get_index_scale_lower_arcsec(index->quads);
        index->indexid = index->quads->indexid;
        index->healpix = index->quads->healpix;
        index->hpnside = index->quads->hpnside;
        index->dimquads = index->quads->dimquads;
        index->nquads = index->quads->numquads;
        index->nstars = index->quads->numstars;
        index->index_jitter = startree_get_jitter(index->starkd);
        if (index->index_jitter == 0.0) {
            index->index_jitter = DEFAULT_INDEX_JITTER;
        }
        index->cutnside = startree_get_cut_nside(index->starkd);
        index->cutnsweep = startree_get_cut_nsweeps(index->starkd);
        index->cutdedup = startree_get_cut_dedup(index->starkd);
        index->cutband = strdup_safe(startree_get_cut_band(index->starkd));
        index->cutmargin = startree_get_cut_margin(index->starkd);
        index_get_missing_cut_params(
            index->indexid,
            index->cutnside == -1 ? &index->cutnside : NULL,
            index->cutnsweep == 0 ? &index->cutnsweep : NULL,
            index->cutdedup == 0 ? &index->cutdedup : NULL,
            index->cutmargin == -1 ? &index->cutmargin : NULL,
            !index->cutband ? &index->cutband : NULL);
        index->circle = qfits_header_getboolean(index->codekd->header, "CIRCLE", 0);
        index->cx_less_than_dx = qfits_header_getboolean(index->codekd->header, "CXDX", FALSE);
        index->meanx_less_than_half = qfits_header_getboolean(index->codekd->header, "CXDXLT1", FALSE);
        pl_append(current->indexes, index);
    }
    if (PyErr_Occurred()) {
        if (current->indexes) {
            pl_remove_all(current->indexes);
            current->indexes = NULL;
        }
        return -1;
    }
    PyObject* logging = PyImport_ImportModule("logging");
    if (!logging) {
        if (current->indexes) {
            pl_remove_all(current->indexes);
            current->indexes = NULL;
        }
        return -1;
    }
    PyObject* message =
        PyUnicode_FromFormat("loaded %zu index file%s", pl_size(current->indexes), pl_size(current->indexes) > 1 ? "s" : "");
    PyObject_CallMethod(logging, "info", "O", message);
    Py_DECREF(message);
    return 0;
}

static PyObject* astrometry_extension_solver_solve(PyObject* self, PyObject* args) {
    PyObject* stars_xs = NULL;
    PyObject* stars_ys = NULL;
    double scale_lower = 0.0;
    double scale_upper = 0.0;
    PyObject* position_hint;
    const char* solve_id;
    int uniformize_index = 0;
    int deduplicate = 0;
    int sip_order = 0;
    int sip_inverse_order = 0;
    int distance_from_quad_bonus = 0;
    double positional_noise_pixels = 0.0;
    double distractor_ratio = 0.0;
    double code_tolerance_l2_distance = 0.0;
    double minimum_quad_size_pixels = 0.0;
    double maximum_quad_size_pixels = 0.0;
    int maximum_quads = 0;
    int maximum_matches = 0;
    int parity = 0;
    PyObject* tune_up_logodds_threshold = NULL;
    double output_logodds_threshold = 0.0;
    PyObject* slices_starts = NULL;
    PyObject* slices_ends = NULL;
    PyObject* logodds_callback = NULL;
    if (!PyArg_ParseTuple(
            args,
            "OOddOsppiiidddddiiiOdOOO",
            &stars_xs,
            &stars_ys,
            &scale_lower,
            &scale_upper,
            &position_hint,
            &solve_id,
            &uniformize_index,
            &deduplicate,
            &sip_order,
            &sip_inverse_order,
            &distance_from_quad_bonus,
            &positional_noise_pixels,
            &distractor_ratio,
            &code_tolerance_l2_distance,
            &minimum_quad_size_pixels,
            &maximum_quad_size_pixels,
            &maximum_quads,
            &maximum_matches,
            &parity,
            &tune_up_logodds_threshold,
            &output_logodds_threshold,
            &slices_starts,
            &slices_ends,
            &logodds_callback)) {
        return NULL;
    }
    astrometry_extension_solver_t* current = (astrometry_extension_solver_t*)self;
    if (scale_lower <= 0.0 || scale_upper <= 0.0 || scale_lower > scale_upper) {
        PyErr_SetString(
            PyExc_TypeError,
            "scale_lower and scale_upper must be strictly positive, and "
            "scale_lower must be smaller than scale_upper");
        return NULL;
    }
    const anbool has_position_hint = position_hint != Py_None;
    double position_hint_ra = 0.0;
    double position_hint_dec = 0.0;
    double position_hint_radius = 0.0;
    if (has_position_hint) {
        if (!PyTuple_Check(position_hint)) {
            PyErr_SetString(PyExc_TypeError, "position_hint must be None or a tuple");
            return NULL;
        }
        if (PyTuple_Size(position_hint) != 3) {
            PyErr_SetString(PyExc_TypeError, "position_hint must have 3 elements");
            return NULL;
        }
        position_hint_ra = PyFloat_AsDouble(PyTuple_GET_ITEM(position_hint, 0));
        if (PyErr_Occurred() || position_hint_ra < 0.0 || position_hint_ra >= 360.0) {
            PyErr_Clear();
            PyErr_SetString(PyExc_TypeError, "position_hint_ra must be a float in the range [0, 360[");
            return NULL;
        }
        position_hint_dec = PyFloat_AsDouble(PyTuple_GET_ITEM(position_hint, 1));
        if (PyErr_Occurred() || position_hint_dec < -90.0 || position_hint_dec > 90.0) {
            PyErr_Clear();
            PyErr_SetString(PyExc_TypeError, "position_hint_dec must be a float in the range [-90, 90]");
            return NULL;
        }
        position_hint_radius = PyFloat_AsDouble(PyTuple_GET_ITEM(position_hint, 2));
        if (PyErr_Occurred() || position_hint_radius < 0) {
            PyErr_Clear();
            PyErr_SetString(PyExc_TypeError, "position_hint_radius must be a float larger than 0");
            return NULL;
        }
    }
    const anbool has_tune_up = tune_up_logodds_threshold != Py_None;
    double tune_up_logodds_threshold_value = 0.0;
    if (has_tune_up) {
        if (!PyFloat_Check(tune_up_logodds_threshold)) {
            PyErr_SetString(PyExc_TypeError, "tune_up_logodds_threshold must be None or a float");
            return NULL;
        }
        tune_up_logodds_threshold_value = PyFloat_AsDouble(tune_up_logodds_threshold);
        if (PyErr_Occurred()) {
            return NULL;
        }
    }
    if (!PyList_Check(stars_xs)) {
        PyErr_SetString(PyExc_TypeError, "stars_xs must be a list");
        return NULL;
    }
    if (!PyList_Check(stars_ys)) {
        PyErr_SetString(PyExc_TypeError, "stars_ys must be a list");
        return NULL;
    }
    const Py_ssize_t stars_size = PyList_GET_SIZE(stars_xs);
    if (stars_size != PyList_GET_SIZE(stars_ys)) {
        PyErr_SetString(PyExc_TypeError, "stars_xs and stars_ys must have the same size");
        return NULL;
    }
    if (stars_size == 0) {
        PyErr_SetString(PyExc_TypeError, "stars_xs cannot be empty");
        return NULL;
    }
    if (!PyList_Check(slices_starts)) {
        PyErr_SetString(PyExc_TypeError, "slices_starts must be a list");
        return NULL;
    }
    if (!PyList_Check(slices_ends)) {
        PyErr_SetString(PyExc_TypeError, "slices_ends must be a list");
        return NULL;
    }
    const Py_ssize_t slices_size = PyList_GET_SIZE(slices_starts);
    if (slices_size != PyList_GET_SIZE(slices_ends)) {
        PyErr_SetString(PyExc_TypeError, "slices_starts and slices_ends must have the same size");
        return NULL;
    }
    if (slices_size == 0) {
        PyErr_SetString(PyExc_TypeError, "slices_starts cannot be empty");
        return NULL;
    }
    int* slices = malloc((size_t)slices_size * 2 * sizeof(double));
    for (Py_ssize_t index = 0; index < slices_size; ++index) {
        slices[index * 2 + 0] = (int)PyLong_AsSize_t(PyList_GET_ITEM(slices_starts, index));
        slices[index * 2 + 1] = (int)PyLong_AsSize_t(PyList_GET_ITEM(slices_ends, index));
        if (PyErr_Occurred() || slices[index * 2 + 0] < 0 || slices[index * 2 + 1] > stars_size
            || slices[index * 2 + 0] >= slices[index * 2 + 1]) {
            free(slices);
            PyErr_SetString(
                PyExc_TypeError,
                "slices_starts and slices_end must contain integers in the range [0, len(stars)[, and slices_starts[i] must be "
                "strictly smaller than slices_ends[i]");
            return NULL;
        }
    }
    PyObject* builtins = PyEval_GetBuiltins();
    if (!builtins) {
        return NULL;
    }
    PyObject* builtins_bool = PyDict_GetItemString(builtins, "bool");
    if (!builtins_bool) {
        return NULL;
    }
    PyObject* sorted_logodds_list = PyList_New(0);
    if (!sorted_logodds_list) {
        return NULL;
    }

    // filter index files using the position and scale hints
    pl* selected_indexes = pl_new(16);
    double maximum_quad_size_arcsec = maximum_quad_size_pixels == 0.0 ? LARGE_VAL : maximum_quad_size_pixels * scale_upper;
    for (size_t position = 0; position < pl_size(current->indexes); ++position) {
        index_t* index = pl_get(current->indexes, position);
        if (index_overlaps_scale_range(index, minimum_quad_size_pixels * scale_lower, maximum_quad_size_arcsec)) {
            if (!has_position_hint
                || (has_position_hint
                    && index_is_within_range(index, position_hint_ra, position_hint_dec, position_hint_radius))) {
                pl_append(selected_indexes, index);
            }
        }
    }

    if (pl_size(selected_indexes) == 0) {
        PyErr_SetString(PyExc_TypeError, "index files do not overlap the provided position and scale hints");
        pl_remove_all(selected_indexes);
        free(slices);
        return NULL;
    }
    PyObject* logging = PyImport_ImportModule("logging");
    if (!logging) {
        pl_remove_all(selected_indexes);
        free(slices);
        return NULL;
    }
    {
        PyObject* message = PyUnicode_FromFormat("solve %s: start", solve_id);
        PyObject_CallMethod(logging, "info", "O", message);
        Py_DECREF(message);
    }
    starxy_t starxy;
    starxy.xlo = LARGE_VAL;
    starxy.xhi = -LARGE_VAL;
    starxy.ylo = LARGE_VAL;
    starxy.yhi = -LARGE_VAL;
    starxy.N = (int)stars_size;
    starxy.x = malloc(sizeof(double) * starxy.N);
    starxy.y = malloc(sizeof(double) * starxy.N);
    starxy.flux = NULL;
    starxy.background = NULL;
    for (Py_ssize_t index = 0; index < stars_size; ++index) {
        starxy.x[index] = PyFloat_AsDouble(PyList_GET_ITEM(stars_xs, index));
        starxy.y[index] = PyFloat_AsDouble(PyList_GET_ITEM(stars_ys, index));
        if (starxy.x[index] < starxy.xlo) {
            starxy.xlo = starxy.x[index];
        }
        if (starxy.x[index] > starxy.xhi) {
            starxy.xhi = starxy.x[index];
        }
        if (starxy.y[index] < starxy.ylo) {
            starxy.ylo = starxy.y[index];
        }
        if (starxy.y[index] > starxy.yhi) {
            starxy.yhi = starxy.y[index];
        }
    }
    if (PyErr_Occurred()) {
        starxy_free_data(&starxy);
        PyErr_Clear();
        PyErr_SetString(PyExc_TypeError, "items in stars_xs and stars_ys must be floats");
        pl_remove_all(selected_indexes);
        Py_DECREF(logging);
        free(slices);
        return NULL;
    }
    callback_context_t context = {
        .solve_id = solve_id,
        .save = PyEval_SaveThread(),
        .builtins_bool = builtins_bool,
        .logging = logging,
        .logodds_callback = logodds_callback,
        .sorted_logodds_list = sorted_logodds_list,
        .solver =
            {
                // input fields
                .indexes = pl_new(1),
                .fieldxy = NULL,
                .pixel_xscale = 0.0,
                .predistort = NULL,
                .fieldxy_orig = &starxy,
                .funits_lower = scale_lower,
                .funits_upper = scale_upper,
                .logratio_toprint =
                    has_tune_up ? MIN(tune_up_logodds_threshold_value, output_logodds_threshold) : output_logodds_threshold,
                .logratio_tokeep = output_logodds_threshold,
                .logratio_totune = has_tune_up ? tune_up_logodds_threshold_value : output_logodds_threshold,
                .record_match_callback = record_match_callback,
                .userdata = NULL,
                .distance_from_quad_bonus = distance_from_quad_bonus,
                .verify_uniformize = uniformize_index,
                .verify_dedup = deduplicate,
                .do_tweak = has_tune_up,
                .tweak_aborder = sip_order,
                .tweak_abporder = sip_inverse_order == 0 ? sip_order : sip_inverse_order,
                .verify_pix = positional_noise_pixels,
                .distractor_ratio = distractor_ratio,
                .codetol = code_tolerance_l2_distance,
                .quadsize_min = minimum_quad_size_pixels,
                .quadsize_max = maximum_quad_size_pixels,
                .startobj = 0,
                .endobj = 0,
                .parity = parity,
                .use_radec = FALSE,
                .centerxyz = {0.0, 0.0, 0.0},
                .r2 = 0.0,
                .logratio_bail_threshold = -1000.0,
                .logratio_stoplooking = LARGE_VAL,
                .maxquads = maximum_quads,
                .maxmatches = maximum_matches,
                .set_crpix = FALSE,
                .set_crpix_center = FALSE,
                .crpix = {0.0, 0.0},
                .mo_template = NULL,
                .timer_callback = timer_callback,

                // callback fields
                .quit_now = FALSE,

                // output fields
                .numtries = 0,
                .nummatches = 0,
                .numscaleok = 0,
                .last_examined_object = 0,
                .num_cxdx_skipped = 0,
                .num_meanx_skipped = 0,
                .num_radec_skipped = 0,
                .num_abscale_skipped = 0,
                .num_verified = 0,

                // internal fields
                .index = NULL,
                .minminAB2 = 0.0,
                .maxmaxAB2 = 0.0,
                .rel_index_noise2 = 0.0,
                .rel_field_noise2 = 0.0,
                .abscale_low = 0.0,
                .abscale_high = 0.0,
                .field_minx = 0.0,
                .field_maxx = 0.0,
                .field_miny = 0.0,
                .field_maxy = 0.0,
                .field_diag = 0.0,
                .cxdx_margin = 0.0,
                .starttime = 0.0,
                .timeused = 0.0,
                .best_logodds = 0.0,
                .best_match = {0},
                .best_index = NULL,
                .best_match_solves = FALSE,
                .have_best_match = FALSE,
                .vf = NULL,
            },
        .output_logodds_threshold = output_logodds_threshold,
        .matches = match_vector_with_capacity(8),
        .error_occured = FALSE,
    };
    {
        err_t* errors_state = errors_get_state();
        errors_state->errfunc = error_callback;
        errors_state->baton = &context;
        errors_state->print_f = NULL;
        errors_state->save = FALSE;
    }
    context.solver.userdata = &context;
    if (has_position_hint) {
        solver_set_radec(&context.solver, position_hint_ra, position_hint_dec, position_hint_radius);
    }
    const size_t selected_indexes_size = pl_size(selected_indexes);
    solver_preprocess_field(&context.solver);
    for (Py_ssize_t slice_index = 0; slice_index < slices_size; ++slice_index) {
        context.solver.startobj = slices[slice_index * 2 + 0];
        context.solver.endobj = slices[slice_index * 2 + 1];
        anbool inner_break = FALSE;
        for (size_t position = 0; position < selected_indexes_size; ++position) {
            index_t* index = pl_get(selected_indexes, position);
            {
                char* index_name = NULL;
                simple_index_name(index, &index_name);
                PyEval_RestoreThread(context.save);
                PyObject* message = PyUnicode_FromFormat(
                    "solve %s: slice=[%d, %d[ (%zd / %zd), index=%s (%zu / %zu)",
                    solve_id,
                    context.solver.startobj,
                    context.solver.endobj,
                    slice_index + 1,
                    slices_size,
                    index_name,
                    position + 1,
                    selected_indexes_size);
                PyObject_CallMethod(logging, "info", "O", message);
                Py_DECREF(message);
                context.save = PyEval_SaveThread();
                free(index_name);
            }
            pl_append(context.solver.indexes, index);
            solver_reset_counters(&context.solver);
            solver_reset_best_match(&context.solver);
            context.solver.have_best_match = TRUE;
            context.solver.best_match.logodds = LARGE_VAL;
            solver_run(&context.solver);
            pl_remove_all(context.solver.indexes);
            context.solver.index = NULL;
            PyEval_RestoreThread(context.save);
            if (PyErr_CheckSignals() != 0) {
                context.error_occured = TRUE;
            }
            if (PyErr_Occurred() || context.error_occured || context.solver.quit_now) {
                context.save = PyEval_SaveThread();
                inner_break = TRUE;
                break;
            }
            context.save = PyEval_SaveThread();
        }
        if (inner_break) {
            break;
        }
    }
    free(slices);
    {
        err_t* errors_state = errors_get_state();
        errors_state->errfunc = NULL;
        errors_state->baton = NULL;
        errors_state->print_f = stderr;
        errors_state->save = FALSE;
    }
    PyEval_RestoreThread(context.save);
    if (PyErr_Occurred() || context.error_occured) {
        if (!PyErr_Occurred()) {
            PyErr_SetNone(PyExc_KeyboardInterrupt);
        }
        context.solver.indexes = NULL;
        context.solver.fieldxy_orig = NULL;
        starxy_free_data(&starxy);
        solver_cleanup(&context.solver);
        pl_remove_all(selected_indexes);
        Py_DECREF(logging);
        Py_DECREF(sorted_logodds_list);
        return NULL;
    }
    PyObject* result = NULL;
    if (context.matches.size > 0) {
        // result = (
        //     stars={
        //         (index_id, star_id): (
        //             ra_deg,
        //             dec_deg,
        //             metadata={key: value, ...}
        //         ),
        //         ...
        //     },
        //     (
        //         (
        //             logodds,
        //             center_ra_deg,
        //             center_dec_deg,
        //             scale,
        //             index_path,
        //             stars_keys=((index_id, star_id), ...),
        //             quad_stars_keys=((index_id, star_id), ...),
        //             wcs_fields={keyword: (value, comment), ...},
        //         ),
        //         ...
        //     ),
        // )
        PyObject* stars = PyDict_New();
        PyObject* matches = PyTuple_New(context.matches.size);
        for (size_t index = 0; index < context.matches.size; ++index) {
            MatchObj* match = context.matches.data + index;
            double ra = 0.0;
            double dec = 0.0;
            xyzarr2radecdeg(match->center, &ra, &dec);
            PyObject* index_id = PyLong_FromLong(match->index->indexid);
            if (match->index->starkd->tagalong == NULL) {
                match->index->starkd->tagalong = get_tagalong(match->index->starkd);
            }
            const int columns =
                match->index->starkd->tagalong == NULL ? 0 : startree_get_tagalong_N_columns(match->index->starkd);
            PyObject* match_stars = PyTuple_New(match->nindex);
            for (int star_id_index = 0; star_id_index < match->nindex; ++star_id_index) {
                PyObject* key = PyTuple_New(2);
                Py_INCREF(index_id);
                PyTuple_SET_ITEM(key, 0, index_id);
                PyTuple_SET_ITEM(key, 1, PyLong_FromLong(match->refstarid[star_id_index]));
                if (PyDict_GetItem(stars, key) == NULL) {
                    PyObject* star = star_to_python_object(
                        match->index->starkd,
                        match->refstarid[star_id_index],
                        match->index->starkd->tagalong != NULL,
                        columns,
                        logging);
                    PyDict_SetItem(stars, key, star);
                }
                PyTuple_SET_ITEM(match_stars, star_id_index, key);
            }
            PyObject* match_quad_stars = PyTuple_New(match->dimquads);
            for (uint8_t quad_index = 0; quad_index < match->dimquads; ++quad_index) {
                kdtree_qres_t* query = kdtree_rangesearch_options(
                    match->index->starkd->tree,
                    match->quadxyz + 3 * (int)quad_index,
                    arcsec2distsq(5.0),
                    KD_OPTIONS_SMALL_RADIUS);
                if (query != NULL && query->nres > 0) {
                    PyObject* key = PyTuple_New(2);
                    Py_INCREF(index_id);
                    PyTuple_SET_ITEM(key, 0, index_id);
                    PyTuple_SET_ITEM(key, 1, PyLong_FromLong(query->inds[0]));
                    PyTuple_SET_ITEM(match_quad_stars, quad_index, key);
                    if (PyDict_GetItem(stars, key) == NULL) {
                        PyObject* star = star_to_python_object(
                            match->index->starkd, query->inds[0], match->index->starkd->tagalong != NULL, columns, logging);
                        PyDict_SetItem(stars, key, star);
                    }
                } else {
                    Py_INCREF(Py_None);
                    PyTuple_SET_ITEM(match_quad_stars, quad_index, Py_None);
                }
                kdtree_free_query(query);
            }
            Py_DECREF(index_id);
            PyObject* wcs_fields = PyDict_New();
            add_wcs_field(wcs_fields, "WCSAXES", PyLong_FromLong(2), "Number of coordinate axes");
            add_wcs_field(wcs_fields, "EQUINOX", PyFloat_FromDouble(2000.0), "Equatorial coordinates definition (yr)");
            add_wcs_field(wcs_fields, "LONPOLE", PyFloat_FromDouble(180.0), "Native longitude of celestial pole (deg)");
            add_wcs_field(wcs_fields, "LATPOLE", PyFloat_FromDouble(0.0), "Native latitude of celestial pole (deg)");
            if (match->sip != NULL && sip_order > 0) {
                add_wcs_field(wcs_fields, "CRVAL1", PyFloat_FromDouble(match->sip->wcstan.crval[0]), "RA of reference point");
                add_wcs_field(wcs_fields, "CRVAL2", PyFloat_FromDouble(match->sip->wcstan.crval[1]), "DEC of reference point");
                add_wcs_field(wcs_fields, "CRPIX1", PyFloat_FromDouble(match->sip->wcstan.crpix[0]), "X reference pixel");
                add_wcs_field(wcs_fields, "CRPIX2", PyFloat_FromDouble(match->sip->wcstan.crpix[1]), "Y reference pixel");
                add_wcs_field(wcs_fields, "CUNIT1", PyUnicode_FromString("deg"), "X pixel scale units");
                add_wcs_field(wcs_fields, "CUNIT2", PyUnicode_FromString("deg"), "Y pixel scale units");
                add_wcs_field(wcs_fields, "CD1_1", PyFloat_FromDouble(match->sip->wcstan.cd[0][0]), "Transformation matrix");
                add_wcs_field(wcs_fields, "CD1_2", PyFloat_FromDouble(match->sip->wcstan.cd[0][1]), "Transformation matrix");
                add_wcs_field(wcs_fields, "CD2_1", PyFloat_FromDouble(match->sip->wcstan.cd[1][0]), "Transformation matrix");
                add_wcs_field(wcs_fields, "CD2_2", PyFloat_FromDouble(match->sip->wcstan.cd[1][1]), "Transformation matrix");
                if (match->sip->wcstan.sin) {
                    add_wcs_field(wcs_fields, "CTYPE1", PyUnicode_FromString("RA---SIN-SIP"), "SIN projection + SIP distortions");
                    add_wcs_field(wcs_fields, "CTYPE2", PyUnicode_FromString("DEC--SIN-SIP"), "SIN projection + SIP distortions");
                } else {
                    add_wcs_field(
                        wcs_fields,
                        "CTYPE1",
                        PyUnicode_FromString("RA---TAN-SIP"),
                        "TAN (gnomonic) projection + SIP distortions");
                    add_wcs_field(
                        wcs_fields,
                        "CTYPE2",
                        PyUnicode_FromString("DEC--TAN-SIP"),
                        "TAN (gnomonic) projection + SIP distortions");
                }
                add_wcs_field(wcs_fields, "A_ORDER", PyLong_FromLong(match->sip->a_order), "Polynomial order, axis 1");
                add_wcs_sip_polynomial(
                    wcs_fields, "A_%i_%i", match->sip->a_order, (double*)match->sip->a, "Polynomial coefficient, axis 1");
                add_wcs_field(wcs_fields, "B_ORDER", PyLong_FromLong(match->sip->b_order), "Polynomial order, axis 2");
                add_wcs_sip_polynomial(
                    wcs_fields, "B_%i_%i", match->sip->b_order, (double*)match->sip->b, "Polynomial coefficient, axis 2");
                add_wcs_field(wcs_fields, "AP_ORDER", PyLong_FromLong(match->sip->ap_order), "Inv polynomial order, axis 1");
                add_wcs_sip_polynomial(
                    wcs_fields, "AP_%i_%i", match->sip->ap_order, (double*)match->sip->ap, "Inv polynomial coefficient, axis 1");
                add_wcs_field(wcs_fields, "BP_ORDER", PyLong_FromLong(match->sip->bp_order), "Inv polynomial order, axis 2");
                add_wcs_sip_polynomial(
                    wcs_fields, "BP_%i_%i", match->sip->bp_order, (double*)match->sip->bp, "Inv polynomial coefficient, axis 2");
            } else {
                add_wcs_field(wcs_fields, "CRVAL1", PyFloat_FromDouble(match->wcstan.crval[0]), "RA of reference point");
                add_wcs_field(wcs_fields, "CRVAL2", PyFloat_FromDouble(match->wcstan.crval[1]), "DEC of reference point");
                add_wcs_field(wcs_fields, "CRPIX1", PyFloat_FromDouble(match->wcstan.crpix[0]), "X reference pixel");
                add_wcs_field(wcs_fields, "CRPIX2", PyFloat_FromDouble(match->wcstan.crpix[1]), "Y reference pixel");
                add_wcs_field(wcs_fields, "CUNIT1", PyUnicode_FromString("deg"), "X pixel scale units");
                add_wcs_field(wcs_fields, "CUNIT2", PyUnicode_FromString("deg"), "Y pixel scale units");
                add_wcs_field(wcs_fields, "CD1_1", PyFloat_FromDouble(match->wcstan.cd[0][0]), "Transformation matrix");
                add_wcs_field(wcs_fields, "CD1_2", PyFloat_FromDouble(match->wcstan.cd[0][1]), "Transformation matrix");
                add_wcs_field(wcs_fields, "CD2_1", PyFloat_FromDouble(match->wcstan.cd[1][0]), "Transformation matrix");
                add_wcs_field(wcs_fields, "CD2_2", PyFloat_FromDouble(match->wcstan.cd[1][1]), "Transformation matrix");
                if (match->wcstan.sin) {
                    add_wcs_field(wcs_fields, "CTYPE1", PyUnicode_FromString("RA---SIN-SIP"), "SIN projection + SIP distortions");
                    add_wcs_field(wcs_fields, "CTYPE2", PyUnicode_FromString("DEC--SIN-SIP"), "SIN projection + SIP distortions");
                } else {
                    add_wcs_field(
                        wcs_fields,
                        "CTYPE1",
                        PyUnicode_FromString("RA---TAN-SIP"),
                        "TAN (gnomonic) projection + SIP distortions");
                    add_wcs_field(
                        wcs_fields,
                        "CTYPE2",
                        PyUnicode_FromString("DEC--TAN-SIP"),
                        "TAN (gnomonic) projection + SIP distortions");
                }
            }
            PyObject* python_match = PyTuple_New(8);
            PyTuple_SET_ITEM(python_match, 0, PyFloat_FromDouble(match->logodds));
            PyTuple_SET_ITEM(python_match, 1, PyFloat_FromDouble(ra));
            PyTuple_SET_ITEM(python_match, 2, PyFloat_FromDouble(dec));
            PyTuple_SET_ITEM(python_match, 3, PyFloat_FromDouble(match->scale));
            PyTuple_SET_ITEM(python_match, 4, PyUnicode_FromString(match->index->indexfn));
            PyTuple_SET_ITEM(python_match, 5, match_stars);
            PyTuple_SET_ITEM(python_match, 6, match_quad_stars);
            PyTuple_SET_ITEM(python_match, 7, wcs_fields);
            PyTuple_SET_ITEM(matches, index, python_match);
            if (has_tune_up) {
                free(match->sip);
            }
            free(match->refradec);
            free(match->fieldxy);
            free(match->fieldxy_orig);
            free(match->theta);
            free(match->matchodds);
            free(match->testperm);
            free(match->refxyz);
            free(match->refxy);
            free(match->refstarid);
            match->sip = NULL;
            match->refradec = NULL;
            match->fieldxy = NULL;
            match->fieldxy_orig = NULL;
            match->tagalong = NULL;
            match->field_tagalong = NULL;
            match->index = NULL;
            match->theta = NULL;
            match->matchodds = NULL;
            match->testperm = NULL;
            match->refxyz = NULL;
            match->refxy = NULL;
            match->refstarid = NULL;
        }
        result = PyTuple_New(2);
        PyTuple_SET_ITEM(result, 0, stars);
        PyTuple_SET_ITEM(result, 1, matches);
    } else {
        Py_INCREF(Py_None);
        result = Py_None;
    }
    match_vector_clear(&context.matches);
    pl_remove_all(context.solver.indexes);
    context.solver.indexes = NULL;
    context.solver.fieldxy_orig = NULL;
    starxy_free_data(&starxy);
    solver_cleanup(&context.solver);
    pl_remove_all(selected_indexes);
    Py_DECREF(logging);
    Py_DECREF(sorted_logodds_list);
    return result;
}

static PyMemberDef astrometry_extension_solver_members[] = {
    {NULL, 0, 0, 0, NULL},
};

static PyTypeObject astrometry_extension_solver_type = {PyVarObject_HEAD_INIT(NULL, 0)};

static PyMethodDef astrometry_extension_solver_methods[] = {
    {"solve", astrometry_extension_solver_solve, METH_VARARGS, NULL},
    {NULL, NULL, 0, NULL},
};

static PyMethodDef astrometry_extension_methods[] = {{NULL, NULL, 0, NULL}};

static struct PyModuleDef astrometry_extension_definition =
    {PyModuleDef_HEAD_INIT, "astrometry_extension", "Astrometry.net core functions wrapper", -1, astrometry_extension_methods};

PyMODINIT_FUNC PyInit_astrometry_extension() {
    PyObject* module = PyModule_Create(&astrometry_extension_definition);
    astrometry_extension_solver_type.tp_name = "astrometry_extension.Solver";
    astrometry_extension_solver_type.tp_basicsize = sizeof(astrometry_extension_solver_t);
    astrometry_extension_solver_type.tp_dealloc = astrometry_extension_solver_dealloc;
    astrometry_extension_solver_type.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE;
    astrometry_extension_solver_type.tp_methods = astrometry_extension_solver_methods;
    astrometry_extension_solver_type.tp_members = astrometry_extension_solver_members;
    astrometry_extension_solver_type.tp_new = astrometry_extension_solver_new;
    astrometry_extension_solver_type.tp_init = astrometry_extension_solver_init;
    PyType_Ready(&astrometry_extension_solver_type);
    PyModule_AddObject(module, "Solver", (PyObject*)&astrometry_extension_solver_type);
    return module;
}
