#include "fitsioutils.h"
#include "mathutil.h"
#include "solver.h"
#include <Python.h>
#include <structmember.h>

#define LARGE_VAL 1e30

typedef struct astrometry_extension_solver_t {
    PyObject_HEAD solver_t* solver;
} astrometry_extension_solver_t;

static void astrometry_extension_solver_dealloc(PyObject* self) {
    astrometry_extension_solver_t* current = (astrometry_extension_solver_t*)self;
    if (current->solver) {
        solver_free(current->solver);
        current->solver = NULL;
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
    current->solver = solver_new();
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
        solver_add_index(current->solver, index);
    }
    if (PyErr_Occurred()) {
        solver_free(current->solver);
        current->solver = NULL;
        return -1;
    }
    PyObject* logging = PyImport_ImportModule("logging");
    if (!logging) {
        solver_free(current->solver);
        current->solver = NULL;
        return -1;
    }
    PyObject* message = PyUnicode_FromFormat(
        "loaded %d index file%s", pl_size(current->solver->indexes), pl_size(current->solver->indexes) > 1 ? "s" : "");
    PyObject_CallMethod(logging, "info", "O", message);
    Py_DECREF(message);
    return 0;
}

typedef struct match_vector_t {
    MatchObj* data;
    size_t size;
    size_t capacity;
} match_vector_t;

static match_vector_t match_vector_with_capacity(size_t capacity) {
    if (capacity == 0) {
        match_vector_t match_vector = {
            NULL,
            0,
            0,
        };
        return match_vector;
    }
    match_vector_t match_vector = {
        malloc(sizeof(MatchObj) * capacity),
        0,
        capacity,
    };
    return match_vector;
}

static void match_vector_clear(match_vector_t* match_vector) {
    if (match_vector->capacity > 0) {
        free(match_vector->data);
        match_vector->data = NULL;
        match_vector->size = 0;
        match_vector->capacity = 0;
    }
}

static void match_vector_reserve(match_vector_t* match_vector, size_t capacity) {
    if (capacity <= match_vector->capacity) {
        return;
    }
    match_vector->data = realloc(match_vector->data, sizeof(MatchObj) * capacity);
    match_vector->capacity = capacity;
}

static void match_vector_push(match_vector_t* match_vector, MatchObj* match) {
    if (match_vector->size == match_vector->capacity) {
        match_vector_reserve(match_vector, match_vector->capacity == 0 ? 1 : match_vector->capacity * 2);
    }
    memcpy(match_vector->data + match_vector->size, match, sizeof(MatchObj));
    match_vector->size += 1;
}

typedef struct callback_context_t {
    const char* solve_id;
    PyThreadState* save;
    PyObject* builtins_bool;
    PyObject* logging;
    PyObject* logodds_callback;
    PyObject* sorted_logodds_list;
    solver_t* solver;
    double output_logodds_threshold;
    match_vector_t matches;
    anbool error_occured;
} callback_context_t;

static anbool record_match_callback(MatchObj* match, void* userdata) {
    callback_context_t* context = (callback_context_t*)userdata;
    verify_hit(
        context->solver->index->starkd,
        context->solver->index->cutnside,
        match,
        match->sip,
        context->solver->vf,
        square(context->solver->verify_pix) + square(context->solver->index->index_jitter / match->scale),
        context->solver->distractor_ratio,
        context->solver->field_maxx,
        context->solver->field_maxy,
        context->solver->logratio_bail_threshold,
        context->solver->logratio_tokeep,
        context->solver->logratio_stoplooking,
        context->solver->distance_from_quad_bonus,
        FALSE);
    if (match->logodds >= context->output_logodds_threshold) {
        double ra = 0.0;
        double dec = 0.0;
        xyzarr2radecdeg(match->center, &ra, &dec);
        char logodds_string[24];
        snprintf(logodds_string, 24, "%g", match->logodds);
        char scale_string[24];
        snprintf(scale_string, 24, "%g", match->scale);
        char ra_string[24];
        snprintf(ra_string, 24, "%g", ra);
        char dec_string[24];
        snprintf(dec_string, 24, "%g", dec);
        match_vector_push(&context->matches, match);
        match->theta = NULL;
        match->matchodds = NULL;
        match->refxyz = NULL;
        match->refxy = NULL;
        match->refstarid = NULL;
        match->testperm = NULL;
        PyEval_RestoreThread(context->save);
        PyObject* message = PyUnicode_FromFormat(
            "solve %s: logodds=%s, matches=%d, conflicts=%d, distractors=%d, "
            "index=%d, ra=%s, dec=%s, scale=%s",
            context->solve_id,
            logodds_string,
            context->matches.data[context->matches.size - 1].nmatch,
            context->matches.data[context->matches.size - 1].nconflict,
            context->matches.data[context->matches.size - 1].ndistractor,
            context->matches.data[context->matches.size - 1].nindex,
            ra_string,
            dec_string,
            scale_string);
        PyObject_CallMethod(context->logging, "info", "O", message);
        Py_DECREF(message);
        anbool inserted_or_error = FALSE;
        for (Py_ssize_t index = 0; index < PyList_Size(context->sorted_logodds_list); ++index) {
            const double value = PyFloat_AsDouble(PyList_GET_ITEM(context->sorted_logodds_list, index));
            if (PyErr_Occurred()) {
                inserted_or_error = TRUE;
                break;
            }
            if (match->logodds > value) {
                PyObject* logodds = PyFloat_FromDouble(match->logodds);
                PyList_Insert(context->sorted_logodds_list, index, logodds);
                Py_DECREF(logodds);
                inserted_or_error = TRUE;
                break;
            }
        }
        if (!inserted_or_error) {
            PyObject* logodds = PyFloat_FromDouble(match->logodds);
            PyList_Append(context->sorted_logodds_list, logodds);
            Py_DECREF(logodds);
        }
        if (PyErr_Occurred()) {
            context->error_occured = TRUE;
            context->solver->quit_now = TRUE;
        } else {
            PyObject* action = PyObject_CallFunction(context->logodds_callback, "O", context->sorted_logodds_list);
            if (PyErr_Occurred()) {
                context->error_occured = TRUE;
                context->solver->quit_now = TRUE;
            } else {
                PyObject* action_as_boolean = PyObject_CallFunction(context->builtins_bool, "O", action);
                if (PyErr_Occurred()) {
                    context->error_occured = TRUE;
                    context->solver->quit_now = TRUE;
                } else if (action_as_boolean == Py_False) {
                    context->solver->quit_now = TRUE;
                }
            }
        }
        const int signal = PyErr_CheckSignals();
        context->save = PyEval_SaveThread();
        if (signal != 0) {
            context->error_occured = TRUE;
            context->solver->quit_now = TRUE;
        }
    } else {
        PyEval_RestoreThread(context->save);
        const int signal = PyErr_CheckSignals();
        context->save = PyEval_SaveThread();
        if (signal != 0) {
            context->error_occured = TRUE;
            context->solver->quit_now = TRUE;
        }
    }
    return FALSE;
}

static time_t timer_callback(void* userdata) {
    callback_context_t* context = (callback_context_t*)userdata;
    PyEval_RestoreThread(context->save);
    const int signal = PyErr_CheckSignals();
    context->save = PyEval_SaveThread();
    if (signal != 0) {
        context->error_occured = TRUE;
        context->solver->quit_now = TRUE;
    }
    return context->solver->quit_now ? 0 : 1;
}

static PyObject*
tagalong_to_python_object(startree_t* tree, int column_index, const char* column_name, int star_id, PyObject* logging) {
    int size = startree_get_tagalong_column_array_size(tree, column_index);
    if (size == 0) {
        Py_INCREF(Py_None);
        return Py_None;
    }
    tfits_type type = startree_get_tagalong_column_fits_type(tree, column_index);
    fitstable_t* table = startree_get_tagalong(tree);
    int row_size = 0;
    void* row = fitstable_read_column_array_inds(table, column_name, type, &star_id, 1, &row_size);
    if (row == NULL) {
        Py_INCREF(Py_None);
        return Py_None;
    }
    PyObject* result = NULL;
    switch (type) {
        case TFITS_BIN_TYPE_D: // double
            if (row_size > 1) {
                result = PyTuple_New(row_size);
                for (int index = 0; index < row_size; ++index) {
                    PyTuple_SET_ITEM(result, index, PyFloat_FromDouble(((double*)row)[index]));
                }
            } else {
                result = PyFloat_FromDouble(*(double*)row);
            }
            break;
        case TFITS_BIN_TYPE_E: // float
            if (row_size > 1) {
                result = PyTuple_New(row_size);
                for (int index = 0; index < row_size; ++index) {
                    PyTuple_SET_ITEM(result, index, PyFloat_FromDouble(((float*)row)[index]));
                }
            } else {
                result = PyFloat_FromDouble(*(float*)row);
            }
            break;
        case TFITS_BIN_TYPE_A: // char
            result = PyUnicode_FromStringAndSize((const char*)row, row_size);
            break;
        case TFITS_BIN_TYPE_I: // i16
            if (row_size > 1) {
                result = PyTuple_New(row_size);
                for (int index = 0; index < row_size; ++index) {
                    PyTuple_SET_ITEM(result, index, PyLong_FromLong(((int16_t*)row)[index]));
                }
            } else {
                result = PyLong_FromLong(*(int16_t*)row);
            }
            break;
        case TFITS_BIN_TYPE_J: // i32
            if (row_size > 1) {
                result = PyTuple_New(row_size);
                for (int index = 0; index < row_size; ++index) {
                    PyTuple_SET_ITEM(result, index, PyLong_FromLong(((int32_t*)row)[index]));
                }
            } else {
                result = PyLong_FromLong(*(int16_t*)row);
            }
            break;
        case TFITS_BIN_TYPE_K: // i64
            if (row_size > 1) {
                result = PyTuple_New(row_size);
                for (int index = 0; index < row_size; ++index) {
                    PyTuple_SET_ITEM(result, index, PyLong_FromLong(((int64_t*)row)[index]));
                }
            } else {
                result = PyLong_FromLong(*(int64_t*)row);
            }
            break;
        default: {
            PyObject* message = PyUnicode_FromFormat("unsupported FITS type %d", type);
            PyObject_CallMethod(logging, "warning", "O", message);
            Py_DECREF(message);
            break;
        }
    }
    free(row);
    if (result == NULL) {
        Py_INCREF(Py_None);
        result = Py_None;
    }
    return result;
}

static PyObject* star_to_python_object(startree_t* tree, int star_id, anbool has_tagalong, int columns, PyObject* logging) {
    double ra = 0.0;
    double dec = 0.0;
    startree_get_radec(tree, star_id, &ra, &dec);
    PyObject* metadata = PyDict_New();
    if (has_tagalong) {
        for (int column_index = 0; column_index < columns; ++column_index) {
            const char* column_name = startree_get_tagalong_column_name(tree, column_index);
            PyObject* value = tagalong_to_python_object(tree, column_index, column_name, star_id, logging);
            PyDict_SetItemString(metadata, column_name, value);
            Py_DECREF(value);
        }
    }
    PyObject* star = PyTuple_New(3);
    PyTuple_SET_ITEM(star, 0, PyFloat_FromDouble(ra));
    PyTuple_SET_ITEM(star, 1, PyFloat_FromDouble(dec));
    PyTuple_SET_ITEM(star, 2, metadata);
    return star;
}

static void add_wcs_field(PyObject* wcs_fields, const char* name, PyObject* value, const char* comment) {
    PyObject* value_and_comment = PyTuple_New(2);
    PyTuple_SET_ITEM(value_and_comment, 0, value);
    PyTuple_SET_ITEM(value_and_comment, 1, PyUnicode_FromString(comment));
    PyDict_SetItemString(wcs_fields, name, value_and_comment);
    Py_DECREF(value_and_comment);
}

static void add_wcs_sip_polynomial(PyObject* wcs_fields, const char* format, int order, const double* data, const char* comment) {
    char name[7];
    // in SIP, data[i * SIP_MAXORDER + j] = 0 if i + j > order
    for (int i = 0; i <= order; ++i) {
        for (int j = 0; i + j <= order; ++j) {
            sprintf(name, format, i, j);
            add_wcs_field(wcs_fields, name, PyFloat_FromDouble(data[i * SIP_MAXORDER + j]), comment);
        }
    }
}

static fitstable_t* get_tagalong(startree_t* star_tree) {
    if (!star_tree->tree->io) {
        return NULL;
    }
    const char* filename = fitsbin_get_filename(star_tree->tree->io);
    if (filename == NULL) {
        return NULL;
    }
    fitstable_t* tag = fitstable_open(filename);
    if (tag == NULL) {
        return NULL;
    }
    int next = fitstable_n_extensions(tag);
    int extension = -1;
    for (int index = 1; index < next; ++index) {
        const qfits_header* header = anqfits_get_header_const(tag->anq, index);
        if (header == NULL) {
            continue;
        }
        char* type = fits_get_dupstring(header, "AN_FILE");
        if (streq(type, AN_FILETYPE_TAGALONG)) {
            free(type);
            extension = index;
            break;
        }
        free(type);
    }
    if (extension == -1) {
        return NULL;
    }
    fitstable_open_extension(tag, extension);
    return tag;
}

static PyObject* astrometry_extension_solver_solve(PyObject* self, PyObject* args) {
    PyObject* stars_xs = NULL;
    PyObject* stars_ys = NULL;
    double scale_lower = 0.0;
    double scale_upper = 0.0;
    PyObject* position_hint;
    const char* solve_id;
    PyObject* tune_up_logodds_threshold = NULL;
    double output_logodds_threshold = 0.0;
    PyObject* logodds_callback = NULL;
    if (!PyArg_ParseTuple(
            args,
            "OOddOsOdO",
            &stars_xs,
            &stars_ys,
            &scale_lower,
            &scale_upper,
            &position_hint,
            &solve_id,
            &tune_up_logodds_threshold,
            &output_logodds_threshold,
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
    PyObject* logging = PyImport_ImportModule("logging");
    if (!logging) {
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
        Py_DECREF(logging);
        return NULL;
    }
    callback_context_t context = {
        .solve_id = solve_id,
        .save = PyEval_SaveThread(),
        .builtins_bool = builtins_bool,
        .logging = logging,
        .logodds_callback = logodds_callback,
        .sorted_logodds_list = sorted_logodds_list,
        .solver = NULL,
        .output_logodds_threshold = output_logodds_threshold,
        .matches = match_vector_with_capacity(8),
        .error_occured = FALSE,
    };
    context.solver = solver_new();
    context.solver->indexes = current->solver->indexes;
    context.solver->fieldxy = NULL;
    context.solver->pixel_xscale = 0.0;
    context.solver->predistort = NULL;
    context.solver->fieldxy_orig = &starxy;
    context.solver->funits_lower = scale_lower;
    context.solver->funits_upper = scale_upper;
    if (has_tune_up) {
        context.solver->logratio_toprint = MIN(tune_up_logodds_threshold_value, output_logodds_threshold);
        context.solver->logratio_tokeep = output_logodds_threshold;
        context.solver->logratio_totune = tune_up_logodds_threshold_value;
        context.solver->do_tweak = TRUE;
    } else {
        context.solver->logratio_toprint = output_logodds_threshold;
        context.solver->logratio_tokeep = output_logodds_threshold;
        context.solver->logratio_totune = output_logodds_threshold;
        context.solver->do_tweak = FALSE;
    }
    context.solver->record_match_callback = record_match_callback;
    context.solver->userdata = &context;
    context.solver->timer_callback = timer_callback;
    if (has_position_hint) {
        solver_set_radec(context.solver, position_hint_ra, position_hint_dec, position_hint_radius);
    }
    // prevent best match copy before running the solver
    context.solver->have_best_match = TRUE;
    context.solver->best_match.logodds = LARGE_VAL;
    solver_run(context.solver);
    PyEval_RestoreThread(context.save);
    if (PyErr_Occurred() || context.error_occured) {
        if (!PyErr_Occurred()) {
            PyErr_SetNone(PyExc_KeyboardInterrupt);
        }
        context.solver->indexes = NULL;
        context.solver->fieldxy_orig = NULL;
        starxy_free_data(&starxy);
        solver_free(context.solver);
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
            if (has_tune_up) {
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
    context.solver->indexes = NULL;
    context.solver->fieldxy_orig = NULL;
    starxy_free_data(&starxy);
    solver_free(context.solver);
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
