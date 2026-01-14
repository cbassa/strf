#include "rftles.h"
#include "satutl.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <csv.h>

// CSV parsing context structure
typedef struct {
    char **fields;           // Current row fields
    size_t field_count;      // Number of fields in current row
    size_t field_capacity;   // Allocated capacity
    size_t current_field;    // Current field being parsed
    
    tle_t *tles;            // Array of TLEs being built
    size_t tle_count;       // Number of TLEs successfully parsed
    size_t tle_capacity;    // Allocated capacity
    
    int header_row;         // Flag: is this the header row?
    int tle0_col;           // Column index for TLE_LINE0
    int tle1_col;           // Column index for TLE_LINE1
    int tle2_col;           // Column index for TLE_LINE2
    int name_col;           // Column index for OBJECT_NAME
    
    size_t row_number;      // Current row number (for error messages)
} csv_parse_ctx_t;

// Callback for CSV field
static void cb_field(void *s, size_t len, void *data) {
    csv_parse_ctx_t *ctx = (csv_parse_ctx_t *)data;
    
    // Ensure we have space for this field
    if (ctx->current_field >= ctx->field_capacity) {
        ctx->field_capacity = ctx->field_capacity ? ctx->field_capacity * 2 : 16;
        ctx->fields = realloc(ctx->fields, ctx->field_capacity * sizeof(char*));
        if (!ctx->fields) {
            fprintf(stderr, "Memory allocation failed for CSV fields\n");
            return;
        }
    }
    
    // Allocate and copy field
    ctx->fields[ctx->current_field] = malloc(len + 1);
    if (ctx->fields[ctx->current_field]) {
        memcpy(ctx->fields[ctx->current_field], s, len);
        ctx->fields[ctx->current_field][len] = '\0';
    }
    
    ctx->current_field++;
    if (ctx->current_field > ctx->field_count) {
        ctx->field_count = ctx->current_field;
    }
}

// Helper function to check if a string looks like a TLE line 1
static int is_tle_line1(const char *str) {
    if (!str || strlen(str) < 2) return 0;
    return (str[0] == '1' && str[1] == ' ');
}

// Helper function to find column index by name
static int find_column(csv_parse_ctx_t *ctx, const char *name) {
    for (size_t i = 0; i < ctx->field_count; i++) {
        if (ctx->fields[i] && strcmp(ctx->fields[i], name) == 0) {
            return (int)i;
        }
    }
    return -1;
}

// Callback for CSV row
static void cb_row(int c, void *data) {
    csv_parse_ctx_t *ctx = (csv_parse_ctx_t *)data;
    ctx->row_number++;
    
    // First row: check if it's a header
    if (ctx->row_number == 1) {
        // Check if this looks like a header by looking for "TLE_LINE1" in the fields
        ctx->tle1_col = find_column(ctx, "TLE_LINE1");
        
        if (ctx->tle1_col >= 0) {
            // Found header columns
            ctx->header_row = 1;
            ctx->tle0_col = find_column(ctx, "TLE_LINE0");
            ctx->tle2_col = find_column(ctx, "TLE_LINE2");
            ctx->name_col = find_column(ctx, "OBJECT_NAME");
            
            if (ctx->tle1_col < 0 || ctx->tle2_col < 0) {
                fprintf(stderr, "Warning: CSV header found but missing required TLE columns\n");
            }
        } else if (ctx->field_count > 39) {
            // No header found, assume fixed positions (columns 37, 38, 39)
            ctx->header_row = 0;
            ctx->tle0_col = 37;
            ctx->tle1_col = 38;
            ctx->tle2_col = 39;
            ctx->name_col = 4;  // OBJECT_NAME typically at column 4
        } else {
            fprintf(stderr, "Warning: CSV format not recognized\n");
            ctx->header_row = 0;
            ctx->tle0_col = -1;
            ctx->tle1_col = -1;
            ctx->tle2_col = -1;
            ctx->name_col = -1;
        }
    }
    
    // Skip header row for data processing
    if (ctx->header_row && ctx->row_number == 1) {
        goto cleanup_row;
    }
    
    // Validate we have the required columns
    if (ctx->tle1_col < 0 || ctx->tle2_col < 0 || 
        (size_t)ctx->tle1_col >= ctx->field_count || 
        (size_t)ctx->tle2_col >= ctx->field_count) {
        goto cleanup_row;
    }
    
    // Get TLE lines
    const char *tle0 = (ctx->tle0_col >= 0 && (size_t)ctx->tle0_col < ctx->field_count) 
                       ? ctx->fields[ctx->tle0_col] : NULL;
    const char *tle1 = ctx->fields[ctx->tle1_col];
    const char *tle2 = ctx->fields[ctx->tle2_col];
    const char *name = (ctx->name_col >= 0 && (size_t)ctx->name_col < ctx->field_count)
                       ? ctx->fields[ctx->name_col] : NULL;
    
    // Validate TLE lines exist and are non-empty
    if (!tle1 || !tle2 || strlen(tle1) < 2 || strlen(tle2) < 2) {
        fprintf(stderr, "Warning: Row %zu: Missing or invalid TLE data", ctx->row_number);
        if (name && strlen(name) > 0) {
            fprintf(stderr, " for satellite '%s'", name);
        }
        fprintf(stderr, "\n");
        goto cleanup_row;
    }
    
    // Basic validation: check if line 1 starts with "1 "
    if (!is_tle_line1(tle1)) {
        fprintf(stderr, "Warning: Row %zu: TLE line 1 doesn't start with '1 '", ctx->row_number);
        if (name && strlen(name) > 0) {
            fprintf(stderr, " for satellite '%s'", name);
        }
        fprintf(stderr, "\n");
        goto cleanup_row;
    }
    
    // Create in-memory file with the three TLE lines
    char tle_buffer[512];
    int offset = 0;
    
    if (tle0 && strlen(tle0) > 0) {
        offset += snprintf(tle_buffer + offset, sizeof(tle_buffer) - offset, "%s\n", tle0);
    }
    offset += snprintf(tle_buffer + offset, sizeof(tle_buffer) - offset, "%s\n", tle1);
    offset += snprintf(tle_buffer + offset, sizeof(tle_buffer) - offset, "%s\n", tle2);
    
    FILE *memfile = fmemopen(tle_buffer, offset, "r");
    if (!memfile) {
        fprintf(stderr, "Warning: Row %zu: Failed to create memory buffer", ctx->row_number);
        if (name && strlen(name) > 0) {
            fprintf(stderr, " for satellite '%s'", name);
        }
        fprintf(stderr, "\n");
        goto cleanup_row;
    }
    
    // Ensure we have space in the TLE array
    if (ctx->tle_count >= ctx->tle_capacity) {
        ctx->tle_capacity = ctx->tle_capacity ? ctx->tle_capacity * 2 : 64;
        ctx->tles = realloc(ctx->tles, ctx->tle_capacity * sizeof(tle_t));
        if (!ctx->tles) {
            fprintf(stderr, "Error: Memory allocation failed for TLE array\n");
            fclose(memfile);
            goto cleanup_row;
        }
    }
    
    // Parse the TLE using existing function
    char satname[ST_SIZE];
    if (read_twoline(memfile, 0, &ctx->tles[ctx->tle_count].orbit, satname) == 0) {
        // Success - store the satellite name
        if (satname[0] != '\0') {
            int satname_len = strlen(satname) + 1;
            ctx->tles[ctx->tle_count].name = malloc(satname_len);
            if (ctx->tles[ctx->tle_count].name) {
                strncpy(ctx->tles[ctx->tle_count].name, satname, satname_len);
            }
        } else {
            ctx->tles[ctx->tle_count].name = NULL;
        }
        ctx->tle_count++;
    } else {
        fprintf(stderr, "Warning: Row %zu: Failed to parse TLE data", ctx->row_number);
        if (name && strlen(name) > 0) {
            fprintf(stderr, " for satellite '%s'", name);
        }
        fprintf(stderr, "\n");
    }
    
    fclose(memfile);
    
cleanup_row:
    // Free all fields for this row
    for (size_t i = 0; i < ctx->field_count; i++) {
        if (ctx->fields[i]) {
            free(ctx->fields[i]);
            ctx->fields[i] = NULL;
        }
    }
    ctx->current_field = 0;
    ctx->field_count = 0;
}

tle_array_t *load_tles_from_csv(char *csvfile) {
    tle_array_t *tle_array;
    
    tle_array = (tle_array_t *)malloc(sizeof(tle_array_t));
    if (tle_array == NULL) {
        return NULL;
    }
    
    tle_array->tles = NULL;
    tle_array->number_of_elements = 0;
    
    // Determine filename
    char filename[1024];
    if (csvfile) {
        strncpy(filename, csvfile, sizeof(filename) - 1);
        filename[sizeof(filename) - 1] = '\0';
    } else {
        char *env = getenv("ST_TLEDIR");
        if (env == NULL || strlen(env) == 0) {
            env = ".";
        }
        snprintf(filename, sizeof(filename), "%s/bulk.csv", env);
    }
    
    // Open CSV file
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "CSV file %s not found\n", filename);
        return tle_array;
    }
    
    // Initialize CSV parser
    struct csv_parser parser;
    if (csv_init(&parser, CSV_STRICT) != 0) {
        fprintf(stderr, "Failed to initialize CSV parser\n");
        fclose(file);
        return tle_array;
    }
    
    // Initialize parsing context
    csv_parse_ctx_t ctx = {0};
    ctx.tle0_col = -1;
    ctx.tle1_col = -1;
    ctx.tle2_col = -1;
    ctx.name_col = -1;
    
    // Parse the CSV file
    char buffer[4096];
    size_t bytes_read;
    
    while ((bytes_read = fread(buffer, 1, sizeof(buffer), file)) > 0) {
        if (csv_parse(&parser, buffer, bytes_read, cb_field, cb_row, &ctx) != bytes_read) {
            fprintf(stderr, "Error parsing CSV: %s\n", csv_strerror(csv_error(&parser)));
            break;
        }
    }
    
    // Finalize parsing
    csv_fini(&parser, cb_field, cb_row, &ctx);
    csv_free(&parser);
    fclose(file);
    
    // Clean up any remaining fields
    if (ctx.fields) {
        for (size_t i = 0; i < ctx.field_capacity; i++) {
            if (ctx.fields[i]) {
                free(ctx.fields[i]);
            }
        }
        free(ctx.fields);
    }
    
    // Transfer parsed TLEs to result
    tle_array->tles = ctx.tles;
    tle_array->number_of_elements = ctx.tle_count;
    
    printf("Loaded %zu orbits from CSV\n", ctx.tle_count);
    
    return tle_array;
}