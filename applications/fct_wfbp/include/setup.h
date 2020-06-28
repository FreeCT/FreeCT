/* FreeCT_wFBP is GPU and CPU CT reconstruction Software */
/* Copyright (C) 2015  John Hoffman */

/* This program is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU General Public License */
/* as published by the Free Software Foundation; either version 2 */
/* of the License, or (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with this program; if not, write to the Free Software */
/* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. */

/* Questions and comments should be directed to */
/* jmhoffman@mednet.ucla.edu with "CTBANGBANG" in the subject line*/

// Non-GPU functions and data types to configure recontructions 
#ifndef setup_h
#define setup_h

#include <recon_structs.h>

// Step 1-3 functions
int configure_paths(struct recon_metadata *mr);
struct recon_params configure_recon_params(char * filename);
struct ct_geom configure_ct_geom(struct recon_metadata *mr);
void configure_reconstruction(struct recon_metadata *mr);
void update_block_info(recon_metadata *mr);
void extract_projections(struct recon_metadata * mr);
void finish_and_cleanup(struct recon_metadata * mr);

#endif
