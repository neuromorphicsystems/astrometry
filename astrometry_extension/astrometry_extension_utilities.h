#include "errors.h"
#include "fitsioutils.h"
#include "mathutil.h"
#include "solver.h"
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <libgen.h>
#include <math.h>
#include <string.h>
#include <structmember.h>

#define LARGE_VAL 1e30

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

static MatchObj* match_vector_pop(match_vector_t* match_vector) {
    if (match_vector->size == 0) {
        return NULL;
    }
    match_vector->size -= 1;
    MatchObj* match = match_vector->data + match_vector->size;
    verify_free_matchobj(match);
    free(match->sip);
    match->sip = NULL;
    return match;
}

static void match_vector_clear(match_vector_t* match_vector) {
    while (match_vector->size > 0) {
        match_vector_pop(match_vector);
    }
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
    if (match_vector->data == NULL) {
        match_vector->data = malloc(sizeof(MatchObj) * capacity);
    } else {
        match_vector->data = realloc(match_vector->data, sizeof(MatchObj) * capacity);
    }
    match_vector->data = realloc(match_vector->data, sizeof(MatchObj) * capacity);
    match_vector->capacity = capacity;
}

static MatchObj* match_vector_push(match_vector_t* match_vector, MatchObj* match) {
    if (match_vector->size == match_vector->capacity) {
        match_vector_reserve(match_vector, match_vector->capacity == 0 ? 1 : match_vector->capacity * 2);
    }
    memcpy(match_vector->data + match_vector->size, match, sizeof(MatchObj));
    MatchObj* inserted_match = match_vector->data + match_vector->size;
    match_vector->size += 1;
    if (match->sip != NULL) {
        sip_t* inserted_sip = malloc(sizeof(sip_t));
        memcpy(inserted_sip, match->sip, sizeof(sip_t));
        inserted_match->sip = inserted_sip;
    }
    inserted_match->refradec = NULL;
    inserted_match->fieldxy = NULL;
    inserted_match->fieldxy_orig = NULL;
    inserted_match->theta = NULL;
    inserted_match->matchodds = NULL;
    inserted_match->testperm = NULL;
    inserted_match->refxyz = NULL;
    inserted_match->refxy = NULL;
    inserted_match->refstarid = NULL;
    return inserted_match;
}

static MatchObj* match_vector_bubble_last(match_vector_t* match_vector) {
    if (match_vector->size == 0) {
        return NULL;
    }
    if (match_vector->size == 1) {
        return match_vector->data;
    }
    float last_logodds = (match_vector->data + match_vector->size - 1)->logodds;
    size_t index = 0;
    for (; index < match_vector->size - 1; ++index) {
        if ((match_vector->data + index)->logodds < last_logodds) {
            MatchObj last;
            memcpy(&last, match_vector->data + match_vector->size - 1, sizeof(MatchObj));
            memmove(
                match_vector->data + index + 1,
                match_vector->data + index,
                (match_vector->size - 1 - index) * sizeof(MatchObj)
            );
            memcpy(match_vector->data + index, &last, sizeof(MatchObj));
            break;
        }
    }
    return match_vector->data + index;
}

static void simple_index_name(index_t* index, char** output) {
    char* index_directory_name = NULL;
    {
        char* match_indexname = strdup(index->indexname);
        char* index_directory = strdup(dirname(match_indexname));
        free(match_indexname);
        match_indexname = NULL;
        index_directory_name = strdup(basename(index_directory));
        free(index_directory);
        index_directory = NULL;
    }
    char* index_name = NULL;
    {
        char* match_indexname = strdup(index->indexname);
        index_name = strdup(basename(match_indexname));
        free(match_indexname);
        match_indexname = NULL;
    }
    const size_t output_size = strlen(index_directory_name) + strlen(index_name) + 4; // +4 for quotes, comma, and null separator
    *output = malloc(output_size);
    sprintf(*output, "\"%s/%s\"", index_directory_name, index_name);
    free(index_name);
    index_name = NULL;
    free(index_directory_name);
    index_directory_name = NULL;
}

typedef struct callback_context_t {
    const char* solve_id;
    PyThreadState* save;
    PyObject* builtins_bool;
    PyObject* logging;
    PyObject* logodds_callback;
    PyObject* logodds_list;
    solver_t solver;
    double output_logodds_threshold;
    match_vector_t matches;
    size_t maximum_matches;
    anbool error_occured;
} callback_context_t;

static anbool record_match_callback(MatchObj* match, void* userdata) {
    callback_context_t* context = (callback_context_t*)userdata;
    MatchObj* inserted_match = match_vector_push(&context->matches, match);
    verify_hit(
        context->solver.index->starkd,
        context->solver.index->cutnside,
        inserted_match,
        inserted_match->sip,
        context->solver.vf,
        square(context->solver.verify_pix) + square(context->solver.index->index_jitter / inserted_match->scale),
        context->solver.distractor_ratio,
        context->solver.field_maxx,
        context->solver.field_maxy,
        context->solver.logratio_bail_threshold,
        context->solver.logratio_tokeep,
        context->solver.logratio_stoplooking,
        context->solver.distance_from_quad_bonus,
        FALSE);
    if (inserted_match->logodds >= context->output_logodds_threshold) {
        inserted_match = match_vector_bubble_last(&context->matches);
        anbool added = TRUE;
        if (context->maximum_matches > 0 && context->matches.size > context->maximum_matches) {
            added = match_vector_pop(&context->matches) != inserted_match;
        }
        if (added) {
            double ra = 0.0;
            double dec = 0.0;
            xyzarr2radecdeg(inserted_match->center, &ra, &dec);
            char logodds_string[24];
            snprintf(logodds_string, 24, "%g", inserted_match->logodds);
            char scale_string[24];
            snprintf(scale_string, 24, "%g", inserted_match->scale);
            char ra_string[24];
            snprintf(ra_string, 24, "%g", ra);
            char dec_string[24];
            snprintf(dec_string, 24, "%g", dec);
            char* index_name = NULL;
            simple_index_name(inserted_match->index, &index_name);
            PyEval_RestoreThread(context->save);
            PyObject* message = PyUnicode_FromFormat(
                "solve %s: logodds=%s, matches=%d, conflicts=%d, distractors=%d, "
                "ra=%s, dec=%s, scale=%s, index=%s",
                context->solve_id,
                logodds_string,
                inserted_match->nmatch,
                inserted_match->nconflict,
                inserted_match->ndistractor,
                ra_string,
                dec_string,
                scale_string,
                index_name);
            PyObject_CallMethod(context->logging, "info", "O", message);
            Py_DECREF(message);
            free(index_name);
            index_name = NULL;
            while (PyList_Size(context->logodds_list) < context->matches.size) {
                PyList_Append(context->logodds_list, Py_None);
            }
            for (Py_ssize_t index = 0; index < PyList_Size(context->logodds_list); ++index) {
                PyList_SetItem(context->logodds_list, index, PyFloat_FromDouble((context->matches.data + index)->logodds));
            }
            if (PyErr_Occurred()) {
                context->error_occured = TRUE;
                context->solver.quit_now = TRUE;
            } else {
                PyObject* action = PyObject_CallFunction(context->logodds_callback, "O", context->logodds_list);
                if (PyErr_Occurred()) {
                    context->error_occured = TRUE;
                    context->solver.quit_now = TRUE;
                } else {
                    PyObject* action_as_boolean = PyObject_CallFunction(context->builtins_bool, "O", action);
                    if (PyErr_Occurred()) {
                        context->error_occured = TRUE;
                        context->solver.quit_now = TRUE;
                    } else if (action_as_boolean == Py_False) {
                        context->solver.quit_now = TRUE;
                    }
                }
            }
            const int signal = PyErr_CheckSignals();
            context->save = PyEval_SaveThread();
            if (signal != 0) {
                context->error_occured = TRUE;
                context->solver.quit_now = TRUE;
            }
        }
    } else {
        match_vector_pop(&context->matches);
        PyEval_RestoreThread(context->save);
        const int signal = PyErr_CheckSignals();
        context->save = PyEval_SaveThread();
        if (signal != 0) {
            context->error_occured = TRUE;
            context->solver.quit_now = TRUE;
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
        context->solver.quit_now = TRUE;
    }
    return context->solver.quit_now ? 0 : 1;
}

static const char* filter_message = "Too few correspondences for the SIP order specified";

static void error_callback(
    void* userdata,
    err_t* error_state,
    const char* module,
    int line,
    const char* func,
    const char* format,
    va_list va) {
    if (strncmp(format, filter_message, strlen(filter_message)) == 0) {
        return;
    }
    callback_context_t* context = (callback_context_t*)userdata;
    PyEval_RestoreThread(context->save);
    PyObject* message = NULL;
    {
        char* new_format = strdup(format);
        const size_t format_size = strlen(format);
        if (format_size > 0 && format[format_size - 1] == '\n') {
            new_format[format_size - 1] = '\0';
        }
        if (line == -1) {
            PyObject* part0 = PyUnicode_FromFormat("%s: ", module);
            PyObject* part1 = PyUnicode_FromFormatV(new_format, va);
            message = PyUnicode_Concat(part0, part1);
            Py_DECREF(part0);
            Py_DECREF(part1);
        } else {
            PyObject* part0 = PyUnicode_FromFormat("%s: ", func);
            PyObject* part1 = PyUnicode_FromFormatV(new_format, va);
            PyObject* part2 = PyUnicode_FromFormat(" (%s, line %d)", module, line);
            PyObject* part01 = PyUnicode_Concat(part0, part1);
            Py_DECREF(part0);
            Py_DECREF(part1);
            message = PyUnicode_Concat(part01, part2);
            Py_DECREF(part01);
            Py_DECREF(part2);
        }
        free(new_format);
        new_format = NULL;
    }
    PyObject_CallMethod(context->logging, "error", "O", message);
    context->save = PyEval_SaveThread();
}

static PyObject* double_to_python_object(double value) {
    if (isnan(value) || isinf(value)) {
        Py_RETURN_NONE;
    }
    return PyFloat_FromDouble(value);
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
                    PyTuple_SET_ITEM(result, index, double_to_python_object(((double*)row)[index]));
                }
            } else {
                result = double_to_python_object(*(double*)row);
            }
            break;
        case TFITS_BIN_TYPE_E: // float
            if (row_size > 1) {
                result = PyTuple_New(row_size);
                for (int index = 0; index < row_size; ++index) {
                    PyTuple_SET_ITEM(result, index, double_to_python_object(((float*)row)[index]));
                }
            } else {
                result = double_to_python_object(*(float*)row);
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
    row = NULL;
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
    char name[9]; // at most 9 characters (AP_10_10 + null terminator)
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
            type = NULL;
            extension = index;
            break;
        }
        free(type);
        type = NULL;
    }
    if (extension == -1) {
        fitstable_close(tag);
        return NULL;
    }
    fitstable_open_extension(tag, extension);
    return tag;
}
